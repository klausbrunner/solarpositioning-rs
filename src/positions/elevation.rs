//! Position-consistent crossings using the geometric event search.

use super::{Model, SolarPositions};
use crate::error::{check_coordinates, check_elevation_angle};
use crate::events::SearchTime;
use crate::{ElevationTarget, Error, Horizon, Location, Result, SolarEvents, grena3, spa};
#[cfg(feature = "chrono")]
use chrono::{DateTime, TimeZone, Utc};

impl SolarPositions {
    /// Finds the next rising solar-centre elevation crossing in `(start, end]`.
    ///
    /// Uses continuous UT Julian dates and the same observation settings as positions.
    /// Requires neither chrono nor allocation. Returns `None` when absent.
    /// Refraction cutoff jumps count as crossings.
    ///
    /// # Errors
    /// Returns an error for invalid coordinates, angle, height, interval or delta T,
    /// times outside the model's UT/TT range, a failed computation, or refraction
    /// too strong to guarantee a monotonic elevation mapping.
    pub fn next_elevation_rise_from_julian(
        &self,
        start: f64,
        end: f64,
        location: Location,
        delta_t: f64,
        target: ElevationTarget,
    ) -> Result<Option<f64>> {
        self.elevation_crossing(start, end, location, delta_t, target, 1.0)
    }

    /// Finds the next setting solar-centre elevation crossing in `(start, end]`.
    ///
    /// # Errors
    /// Has the same errors as [`Self::next_elevation_rise_from_julian`].
    pub fn next_elevation_set_from_julian(
        &self,
        start: f64,
        end: f64,
        location: Location,
        delta_t: f64,
        target: ElevationTarget,
    ) -> Result<Option<f64>> {
        self.elevation_crossing(start, end, location, delta_t, target, -1.0)
    }

    /// Finds the next rising solar-centre elevation crossing in `(start, end]`.
    ///
    /// Returns the start's time zone. Cutoff jumps count as crossings.
    ///
    /// # Errors
    /// Has the same errors as [`Self::next_elevation_rise_from_julian`].
    #[cfg(feature = "chrono")]
    #[cfg_attr(docsrs, doc(cfg(feature = "chrono")))]
    pub fn next_elevation_rise<Tz: TimeZone>(
        &self,
        start: &DateTime<Tz>,
        end: &DateTime<Tz>,
        location: Location,
        delta_t: f64,
        target: ElevationTarget,
    ) -> Result<Option<DateTime<Tz>>> {
        let zone = start.timezone();
        self.elevation_crossing(
            start.with_timezone(&Utc),
            end.with_timezone(&Utc),
            location,
            delta_t,
            target,
            1.0,
        )
        .map(|event| event.map(|time| time.with_timezone(&zone)))
    }

    /// Finds the next setting solar-centre elevation crossing in `(start, end]`.
    ///
    /// # Errors
    /// Has the same errors as [`Self::next_elevation_rise_from_julian`].
    #[cfg(feature = "chrono")]
    #[cfg_attr(docsrs, doc(cfg(feature = "chrono")))]
    pub fn next_elevation_set<Tz: TimeZone>(
        &self,
        start: &DateTime<Tz>,
        end: &DateTime<Tz>,
        location: Location,
        delta_t: f64,
        target: ElevationTarget,
    ) -> Result<Option<DateTime<Tz>>> {
        let zone = start.timezone();
        self.elevation_crossing(
            start.with_timezone(&Utc),
            end.with_timezone(&Utc),
            location,
            delta_t,
            target,
            -1.0,
        )
        .map(|event| event.map(|time| time.with_timezone(&zone)))
    }

    fn elevation_crossing<T: SearchTime>(
        self,
        start: T,
        end: T,
        location: Location,
        delta_t: f64,
        target: ElevationTarget,
        direction: f64,
    ) -> Result<Option<T>> {
        check_coordinates(location.latitude, location.longitude)?;
        check_elevation_angle(target.angle)?;
        if !target.height.is_finite()
            || (matches!(self.model, Model::Grena3) && target.height != 0.0)
        {
            return Err(Error::InvalidHeight {
                value: target.height,
            });
        }
        let events = SolarEvents::with_range(
            |time, location| match self.model {
                Model::Spa => spa::event_position_at_height(time, location, target.height),
                Model::Grena3 => grena3::event_position(time, location),
            },
            self.model.range(),
        );
        events.check_interval(start, end, delta_t)?;
        let Some(angle) = self.geometric_threshold(target)? else {
            return Ok(None);
        };
        events.crossing(
            start,
            end,
            location,
            delta_t,
            Horizon::Custom(angle),
            direction,
        )
    }

    fn geometric_threshold(self, target: ElevationTarget) -> Result<Option<f64>> {
        let Some(correction) = target.refraction else {
            return Ok(Some(target.angle));
        };
        let (cutoff, slope_bound) = match self.model {
            Model::Spa => (
                spa::REFRACTION_CUTOFF,
                spa::refraction_slope_bound(correction),
            ),
            Model::Grena3 => (
                grena3::REFRACTION_CUTOFF,
                grena3::refraction_slope_bound(correction),
            ),
        };
        if slope_bound >= 1.0 {
            return Err(Error::ComputationError {
                message: "refraction too strong to guarantee monotonic elevation crossings",
            });
        }
        if target.angle <= cutoff {
            return Ok(Some(target.angle));
        }
        let apparent = |angle: f64| {
            angle
                + match self.model {
                    Model::Spa => spa::refraction_correction(angle, correction),
                    Model::Grena3 => {
                        grena3::refraction_correction(angle.to_radians(), correction).to_degrees()
                    }
                }
        };
        if target.angle <= apparent(cutoff) {
            return Ok(Some(cutoff));
        }
        if target.angle > apparent(90.0) {
            return Ok(None);
        }
        let (mut low, mut high) = (cutoff, 90.0);
        // Invert only the continuous active branch, once per search.
        for _ in 0..48 {
            let middle = low + (high - low) / 2.0;
            if apparent(middle) < target.angle {
                low = middle;
            } else {
                high = middle;
            }
        }
        Ok(Some(high))
    }
}
