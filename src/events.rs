//! Solar events found from precise positions, with SPA as the default model.
//!
//! Bounded searches exclude their start and include their end. Pass a returned time
//! as the next start to continue without repeating an event. Rise, set and upper
//! meridian transit are independent. A tangency alone is not a crossing.
//!
//! Positions are unrefracted, topocentric solar-centre positions at sea level.
//! [`Horizon::SunriseSunset`] includes the conventional allowance for refraction and
//! solar radius; other horizons use their specified geometric elevation. Brackets
//! are refined to one millisecond, which does not imply equivalent observed-event
//! accuracy. Crossings closer together need not be distinguished.
//!
//! Numeric times are continuous UT Julian dates. Chrono helpers use a proleptic
//! Gregorian calendar, approximating UT1 with UTC. Delta T (TT minus UT1, in seconds)
//! is constant during a search. Both UT and TT must remain in the provider's range.

use crate::error::{check_coordinates, check_elevation_angle};
use crate::math::{abs, copysign, cos, sin};
use crate::time::JulianDate;
use crate::{Error, Horizon, Location, Result, grena3, spa};
use core::ops::RangeInclusive;

#[cfg(feature = "alloc")]
use {alloc::vec::Vec, core::ops::Range};
#[cfg(feature = "chrono")]
mod chrono;

const TIME_TOLERANCE: f64 = 1.0 / 3_600_000.0;
// Conservative allowance for Julian-date quantisation and floating-point evaluation.
const ROUNDING_ERROR: f64 = 1e-8;

/// Coordinates needed by the event search, in degrees.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EventPosition {
    /// Unrefracted topocentric solar-centre elevation at sea level.
    pub elevation: f64,
    /// Geocentric local hour angle, positive west; any full-turn wrapping is accepted.
    pub hour_angle: f64,
}

/// Function-pointer provider used by the built-in calculators.
///
/// Custom calculators can also accept closures through [`SolarEvents::with_provider`].
pub type PositionFn = fn(JulianDate, Location) -> Result<EventPosition>;

/// A reusable solar-event calculator, defaulting to SPA.
///
/// This type owns its provider. Custom providers can be functions or closures; a
/// calculator is `Send` and `Sync` when its provider is.
#[derive(Debug, Clone, Copy)]
pub struct SolarEvents<P = PositionFn> {
    provider: P,
    min_time: f64,
    max_time: f64,
}

impl SolarEvents {
    /// Creates a calculator using SPA, for years -2000 through 6000.
    #[must_use]
    pub fn new() -> Self {
        Self {
            provider: spa::event_position,
            min_time: year_start(-2000),
            max_time: year_start(6001),
        }
    }

    /// Creates a calculator using Grena3, for years 2010 through 2110.
    #[must_use]
    pub fn grena3() -> Self {
        Self {
            provider: grena3::event_position,
            min_time: year_start(2010),
            max_time: year_start(2111),
        }
    }
}

impl Default for SolarEvents {
    fn default() -> Self {
        Self::new()
    }
}

impl<P> SolarEvents<P>
where
    P: Fn(JulianDate, Location) -> Result<EventPosition>,
{
    /// Creates a calculator with a custom solar-position model and inclusive year range.
    ///
    /// The provider receives continuous UT/TT Julian dates and a [`Location`]
    /// in degrees (longitude positive east). It must return finite, deterministic
    /// coordinates. This is a solar-model extension, not an arbitrary-object search:
    /// the sines of elevation and hour angle must be continuous, with absolute
    /// second derivatives in hours below `0.1 * abs(cos(latitude)) + 0.0001` and
    /// `0.1`, respectively. Hour angle is independent of latitude and increases
    /// westwards through zero at upper transit.
    ///
    /// # Errors
    /// Returns an error for an empty year range or years outside -2000 through 6000.
    pub fn with_provider(provider: P, years: RangeInclusive<i32>) -> Result<Self> {
        if years.is_empty() || *years.start() < -2000 || *years.end() > 6000 {
            return Err(invalid_time("invalid provider year range"));
        }
        Ok(Self {
            provider,
            min_time: year_start(*years.start()),
            max_time: year_start(*years.end() + 1),
        })
    }

    /// Finds the next rising crossing in `(start, end]`, using UT Julian dates.
    ///
    /// Returns `None` when absent. This method needs neither chrono nor allocation.
    ///
    /// # Errors
    /// Returns an error for invalid coordinates, horizon, interval or delta T, times
    /// outside the provider's UT/TT range, or a failed/non-finite provider result.
    pub fn next_rise_from_julian(
        &self,
        start: f64,
        end: f64,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Option<f64>> {
        self.crossing(start, end, location, delta_t, horizon, 1.0)
    }

    /// Finds the next setting crossing in `(start, end]`, using UT Julian dates.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_rise_from_julian`].
    pub fn next_set_from_julian(
        &self,
        start: f64,
        end: f64,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Option<f64>> {
        self.crossing(start, end, location, delta_t, horizon, -1.0)
    }

    /// Finds the next upper meridian transit in `(start, end]`, using UT Julian dates.
    ///
    /// Transit is independent of latitude and horizon, and need not coincide with
    /// maximum elevation.
    ///
    /// # Errors
    /// Returns an error for invalid longitude, interval or delta T, times outside
    /// the provider's UT/TT range, or a failed/non-finite provider result.
    pub fn next_transit_from_julian(
        &self,
        start: f64,
        end: f64,
        longitude: f64,
        delta_t: f64,
    ) -> Result<Option<f64>> {
        self.transit(start, end, longitude, delta_t)
    }

    /// Collects all events in `[start, end)`, using UT Julian dates.
    ///
    /// Use consecutive UTC midnights for a UTC calendar day, or supply other date
    /// boundaries. Unlike bounded next-event searches, collection includes its
    /// start and excludes its end (to one-millisecond resolution).
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_rise_from_julian`].
    #[cfg(feature = "alloc")]
    #[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
    pub fn for_interval_from_julian(
        &self,
        start: f64,
        end: f64,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Events<f64>> {
        let interval = start..end;
        let (initial, transits) = self.prepare_interval(&interval, location, delta_t)?;
        self.horizon_events(&interval, location, delta_t, horizon, initial, transits)
    }

    fn position(&self, jd: f64, location: Location, delta_t: f64) -> Result<EventPosition> {
        let position = (self.provider)(JulianDate::new(jd, delta_t)?, location)?;
        if !position.elevation.is_finite() || !position.hour_angle.is_finite() {
            return Err(Error::ComputationError {
                message: "non-finite solar position",
            });
        }
        Ok(position)
    }

    #[allow(clippy::suboptimal_flops)]
    fn crossing<T: SearchTime>(
        &self,
        start: T,
        end: T,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
        direction: f64,
    ) -> Result<Option<T>> {
        check_coordinates(location.latitude, location.longitude)?;
        check_elevation_angle(horizon.elevation_angle())?;
        let horizon = sin(horizon.elevation_angle().to_radians());
        // Daily rotation contributes about (2*pi/24)^2*cos(latitude) = 0.069*cos(latitude).
        // Round up to 0.1, with a floor for slower solar motion and parallax at the poles.
        let curvature = 0.1 * abs(cos(location.latitude.to_radians())) + 0.0001;
        self.search(start, end, delta_t, curvature, |jd| {
            Ok(direction
                * (sin(self.position(jd, location, delta_t)?.elevation.to_radians()) - horizon))
        })
    }

    fn transit<T: SearchTime>(
        &self,
        start: T,
        end: T,
        longitude: f64,
        delta_t: f64,
    ) -> Result<Option<T>> {
        check_coordinates(0.0, longitude)?;
        self.search(start, end, delta_t, 0.1, |jd| {
            Ok(sin(self
                .position(
                    jd,
                    Location {
                        latitude: 0.0,
                        longitude,
                    },
                    delta_t,
                )?
                .hour_angle
                .to_radians()))
        })
    }

    fn search<T: SearchTime>(
        &self,
        start: T,
        end: T,
        delta_t: f64,
        curvature: f64,
        position: impl Fn(f64) -> Result<f64>,
    ) -> Result<Option<T>> {
        self.check_interval(start, end, delta_t)?;
        let hours = start.hours_until(end);
        if hours == 0.0 {
            return Ok(None);
        }
        // Evaluate exactly the time that will be returned, including timestamp rounding.
        // Otherwise resuming at a returned crossing can rediscover it.
        let search = CrossingSearch {
            function: |t| position(start.at(end, t, hours).julian_day()),
            curvature,
        };
        let crossing = search.find(0.0, hours, search.value(0.0)?, search.value(hours)?)?;
        Ok(crossing.map(|t| start.at(end, t, hours)))
    }

    fn check_interval<T: SearchTime>(&self, start: T, end: T, delta_t: f64) -> Result<()> {
        if !delta_t.is_finite() || end < start {
            return Err(invalid_time("invalid search interval or delta_t"));
        }
        for time in [start, end] {
            let tt = time.julian_day() + delta_t / 86400.0;
            if !time.in_range(self.min_time, self.max_time)
                || !(self.min_time..=self.max_time).contains(&tt)
            {
                return Err(invalid_time(
                    "search interval outside provider's supported years",
                ));
            }
        }
        Ok(())
    }

    #[cfg(feature = "alloc")]
    fn collect<T: SearchTime>(
        &self,
        start: T,
        end: T,
        delta_t: f64,
        next: impl Fn(T) -> Result<Option<T>>,
    ) -> Result<Vec<T>> {
        // Include boundary crossings, without looking outside the provider's UT/TT range.
        let lookback = start.shift_hours(-TIME_TOLERANCE);
        let mut cursor = if self.check_interval(lookback, end, delta_t).is_ok() {
            lookback
        } else {
            start
        };
        let mut events = Vec::new();
        while let Some(time) = next(cursor)? {
            if time >= end {
                break;
            }
            if time >= start {
                events.push(time);
            }
            cursor = time;
        }
        Ok(events)
    }

    // Shared work for all horizons: initial sine of elevation and meridian transits.
    #[cfg(feature = "alloc")]
    fn prepare_interval<T: SearchTime>(
        &self,
        interval: &Range<T>,
        location: Location,
        delta_t: f64,
    ) -> Result<(f64, Vec<T>)> {
        check_coordinates(location.latitude, location.longitude)?;
        self.check_interval(interval.start, interval.end, delta_t)?;
        let transits = self.collect(interval.start, interval.end, delta_t, |t| {
            self.transit(t, interval.end, location.longitude, delta_t)
        })?;
        let initial = sin(self
            .position(interval.start.julian_day(), location, delta_t)?
            .elevation
            .to_radians());
        Ok((initial, transits))
    }

    #[cfg(feature = "alloc")]
    fn horizon_events<T: SearchTime>(
        &self,
        interval: &Range<T>,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
        initial: f64,
        transits: Vec<T>,
    ) -> Result<Events<T>> {
        check_elevation_angle(horizon.elevation_angle())?;
        let difference = initial - sin(horizon.elevation_angle().to_radians());
        let state_at_start = if abs(difference) <= ROUNDING_ERROR {
            HorizonState::OnHorizon
        } else if difference > 0.0 {
            HorizonState::Above
        } else {
            HorizonState::Below
        };
        Ok(Events {
            start: interval.start,
            end: interval.end,
            state_at_start,
            rises: self.collect(interval.start, interval.end, delta_t, |t| {
                self.crossing(t, interval.end, location, delta_t, horizon, 1.0)
            })?,
            transits,
            sets: self.collect(interval.start, interval.end, delta_t, |t| {
                self.crossing(t, interval.end, location, delta_t, horizon, -1.0)
            })?,
        })
    }
}

/// Position relative to a horizon at the beginning of an interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HorizonState {
    /// Above the selected horizon.
    Above,
    /// Below the selected horizon.
    Below,
    /// Within the numerical rounding allowance of the horizon.
    OnHorizon,
}

/// All solar events in an interval with inclusive start and exclusive end.
///
/// Each list is chronological and can be empty or contain multiple events.
/// A skipped calendar date has equal start and end and empty lists.
#[cfg(feature = "alloc")]
#[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Events<T> {
    /// Inclusive interval start.
    pub start: T,
    /// Exclusive interval end.
    pub end: T,
    /// Solar-centre state relative to the horizon at start.
    pub state_at_start: HorizonState,
    /// Rising crossings.
    pub rises: Vec<T>,
    /// Upper meridian transits, independent of horizon.
    pub transits: Vec<T>,
    /// Setting crossings.
    pub sets: Vec<T>,
}

#[cfg(feature = "alloc")]
impl<T: PartialOrd> Events<T> {
    /// No horizon crossings and a start above the horizon; false for an empty interval.
    #[must_use]
    pub fn always_above(&self) -> bool {
        self.end > self.start
            && self.state_at_start == HorizonState::Above
            && self.rises.is_empty()
            && self.sets.is_empty()
    }

    /// No horizon crossings and a start below the horizon; false for an empty interval.
    #[must_use]
    pub fn always_below(&self) -> bool {
        self.end > self.start
            && self.state_at_start == HorizonState::Below
            && self.rises.is_empty()
            && self.sets.is_empty()
    }
}

const fn invalid_time(message: &'static str) -> Error {
    Error::InvalidDateTime { message }
}

// Proleptic Gregorian January 1, including astronomical year zero.
fn year_start(year: i32) -> f64 {
    let y = year - 1;
    1_721_425.5 + f64::from(365 * y + y.div_euclid(4) - y.div_euclid(100) + y.div_euclid(400))
}

// Keep both numeric and chrono searches on the same algorithm, while evaluating
// the exact representation returned to their respective callers.
trait SearchTime: Copy + PartialOrd {
    fn julian_day(self) -> f64;
    fn hours_until(self, end: Self) -> f64;
    fn shift_hours(self, hours: f64) -> Self;
    fn in_range(self, min: f64, max: f64) -> bool;

    fn at(self, end: Self, hours: f64, duration_hours: f64) -> Self {
        // Preserve the inclusive endpoint exactly through the floating-point round trip.
        if hours == duration_hours {
            return end;
        }
        let time = self.shift_hours(hours);
        if time > end { end } else { time }
    }
}

impl SearchTime for f64 {
    fn julian_day(self) -> f64 {
        self
    }
    fn hours_until(self, end: Self) -> f64 {
        (end - self) * 24.0
    }
    fn shift_hours(self, hours: f64) -> Self {
        self + hours / 24.0
    }
    fn in_range(self, min: f64, max: f64) -> bool {
        (min..=max).contains(&self)
    }
}

struct CrossingSearch<F> {
    function: F,
    curvature: f64,
}

impl<F: Fn(f64) -> Result<f64>> CrossingSearch<F> {
    fn value(&self, time: f64) -> Result<f64> {
        (self.function)(time)
    }

    fn find(&self, start: f64, end: f64, a: f64, b: f64) -> Result<Option<f64>> {
        let width = end - start;
        let slope = (b - a) / width;
        let slope_error = self.curvature * width / 2.0 + 2.0 * ROUNDING_ERROR / width;
        // Standard interpolation bounds: slope error M*w/2, chord error M*w*w/8.
        // https://dlmf.nist.gov/3.3.E5. Discard only intervals that cannot rise
        // through zero; otherwise subdivide, searching the earlier half first.
        if slope + slope_error <= 0.0
            || ((a > 0.0) == (b > 0.0)
                && abs(a).min(abs(b)) > self.curvature * width * width / 8.0 + ROUNDING_ERROR)
        {
            return Ok(None);
        }
        if a < 0.0 && b >= 0.0 && slope > slope_error {
            return self.refine(start, end, a, b).map(Some);
        }
        if width <= TIME_TOLERANCE {
            return Ok((a < 0.0 && b > 0.0).then_some(end));
        }
        let middle = start + width / 2.0;
        let m = self.value(middle)?;
        if let Some(left) = self.find(start, middle, a, m)? {
            return Ok(Some(left));
        }
        self.find(middle, end, m, b)
    }

    #[allow(clippy::while_float)] // The projection bounds the number of refinement steps.
    fn refine(&self, mut start: f64, mut end: f64, mut a: f64, mut b: f64) -> Result<f64> {
        // ITP (Oliveira & Takahashi): https://doi.org/10.1145/3423597.
        // Use k1 = 0.2 / initial width, k2 = 2, n0 = 1. The projection allows
        // at most one extra iteration over bisection while favouring interpolation.
        let scale = 0.2 / (end - start);
        // Maximum remaining width after a step; halve the allowance each time.
        let mut max_width = TIME_TOLERANCE;
        while max_width < end - start {
            max_width *= 2.0;
        }
        while end - start > TIME_TOLERANCE {
            let width = end - start;
            let middle = start + width / 2.0;
            // Equal values can occur on a rounded zero; use the midpoint then.
            let interpolated = if a == b {
                middle
            } else {
                start - a * width / (b - a)
            };
            let towards_middle = middle - interpolated;
            let truncated = interpolated
                + copysign(
                    (scale * width * width).min(abs(towards_middle)),
                    towards_middle,
                );
            let radius = (max_width - width / 2.0).max(0.0);
            let trial = middle + (truncated - middle).clamp(-radius, radius);
            let value = self.value(trial)?;
            if value > 0.0 {
                end = trial;
                b = value;
            } else {
                start = trial;
                a = value;
            }
            max_width /= 2.0;
        }
        // Keep the later endpoint, so the next search cannot rediscover this crossing.
        Ok(end)
    }
}
