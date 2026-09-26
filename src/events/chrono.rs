#[cfg(feature = "alloc")]
use super::invalid_time;
use super::{EventPosition, SearchTime, SolarEvents};
use crate::{
    Horizon, Location, Result,
    time::{JulianDate, datetime_to_julian},
};
use ::chrono::{DateTime, Duration, TimeZone, Utc};
#[cfg(feature = "alloc")]
use {super::Events, ::chrono::NaiveDate, alloc::vec::Vec, core::ops::Range};

impl<P> SolarEvents<P>
where
    P: Fn(JulianDate, Location) -> Result<EventPosition>,
{
    /// Finds the next rising crossing in `(start, end]`, returning the start's time zone.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_rise_from_julian`].
    pub fn next_rise<Tz: TimeZone>(
        &self,
        start: &DateTime<Tz>,
        end: &DateTime<Tz>,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Option<DateTime<Tz>>> {
        let zone = start.timezone();
        self.crossing(
            start.with_timezone(&Utc),
            end.with_timezone(&Utc),
            location,
            delta_t,
            horizon,
            1.0,
        )
        .map(|event| event.map(|time| time.with_timezone(&zone)))
    }

    /// Finds the next setting crossing in `(start, end]`, returning the start's time zone.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_rise_from_julian`].
    pub fn next_set<Tz: TimeZone>(
        &self,
        start: &DateTime<Tz>,
        end: &DateTime<Tz>,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Option<DateTime<Tz>>> {
        let zone = start.timezone();
        self.crossing(
            start.with_timezone(&Utc),
            end.with_timezone(&Utc),
            location,
            delta_t,
            horizon,
            -1.0,
        )
        .map(|event| event.map(|time| time.with_timezone(&zone)))
    }

    /// Finds the next upper meridian transit in `(start, end]`.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_transit_from_julian`].
    pub fn next_transit<Tz: TimeZone>(
        &self,
        start: &DateTime<Tz>,
        end: &DateTime<Tz>,
        longitude: f64,
        delta_t: f64,
    ) -> Result<Option<DateTime<Tz>>> {
        let zone = start.timezone();
        self.transit(
            start.with_timezone(&Utc),
            end.with_timezone(&Utc),
            longitude,
            delta_t,
        )
        .map(|event| event.map(|time| time.with_timezone(&zone)))
    }

    /// Collects all events in a local calendar date, including its start and excluding its end.
    ///
    /// Handles daylight-saving changes, repeated dates and skipped dates. Returned
    /// times use `zone`. Pass [`Horizon::SunriseSunset`] for standard rise/set times.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_rise_from_julian`].
    /// Also returns an error if the time zone supplies no valid local time within
    /// a day of a date boundary, or the following calendar date is out of range.
    #[cfg(feature = "alloc")]
    #[cfg_attr(docsrs, doc(cfg(all(feature = "chrono", feature = "alloc"))))]
    pub fn for_date<Tz: TimeZone>(
        &self,
        date: NaiveDate,
        zone: &Tz,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Events<DateTime<Tz>>> {
        let interval = date_interval(date, zone)?;
        let (initial, transits) = self.prepare_interval(&interval, location, delta_t)?;
        self.horizon_events(&interval, location, delta_t, horizon, initial, transits)
            .map(|events| events.into_timezone(zone))
    }

    /// Collects a date's events for several horizons, sharing the transit calculation.
    ///
    /// Results follow input order, including repeated horizons. Each horizon has its
    /// own initial state and crossing lists.
    ///
    /// # Errors
    /// Has the same errors as [`Self::for_date`].
    #[cfg(feature = "alloc")]
    #[cfg_attr(docsrs, doc(cfg(all(feature = "chrono", feature = "alloc"))))]
    pub fn for_date_multiple<Tz: TimeZone>(
        &self,
        date: NaiveDate,
        zone: &Tz,
        location: Location,
        delta_t: f64,
        horizons: impl IntoIterator<Item = Horizon>,
    ) -> Result<Vec<(Horizon, Events<DateTime<Tz>>)>> {
        let interval = date_interval(date, zone)?;
        let (initial, transits) = self.prepare_interval(&interval, location, delta_t)?;
        horizons
            .into_iter()
            .map(|horizon| {
                self.horizon_events(
                    &interval,
                    location,
                    delta_t,
                    horizon,
                    initial,
                    transits.clone(),
                )
                .map(|events| (horizon, events.into_timezone(zone)))
            })
            .collect()
    }
}

#[cfg(feature = "alloc")]
impl Events<DateTime<Utc>> {
    fn into_timezone<Tz: TimeZone>(self, zone: &Tz) -> Events<DateTime<Tz>> {
        let convert = |time: DateTime<Utc>| time.with_timezone(zone);
        Events {
            start: convert(self.start),
            end: convert(self.end),
            state_at_start: self.state_at_start,
            rises: self.rises.into_iter().map(convert).collect(),
            transits: self.transits.into_iter().map(convert).collect(),
            sets: self.sets.into_iter().map(convert).collect(),
        }
    }
}

#[cfg(feature = "alloc")]
fn date_interval<Tz: TimeZone>(date: NaiveDate, zone: &Tz) -> Result<Range<DateTime<Utc>>> {
    let next_date = date
        .succ_opt()
        .ok_or_else(|| invalid_time("calendar date out of range"))?;
    let start = start_of_date(date, zone)?.with_timezone(&Utc);
    let end = start_of_date(next_date, zone)?.with_timezone(&Utc);
    Ok(start..end)
}

impl SearchTime for DateTime<Utc> {
    fn julian_day(self) -> f64 {
        datetime_to_julian(&self)
    }

    fn hours_until(self, end: Self) -> f64 {
        let duration = end - self;
        (duration.num_seconds() as f64 + f64::from(duration.subsec_nanos()) / 1e9) / 3600.0
    }

    fn shift_hours(self, hours: f64) -> Self {
        let seconds = hours * 3600.0;
        let whole_seconds = seconds as i64;
        self + Duration::seconds(whole_seconds)
            + Duration::nanoseconds(((seconds - whole_seconds as f64) * 1e9) as i64)
    }

    fn in_range(self, min: f64, max: f64) -> bool {
        let min_seconds = ((min - 2_440_587.5) * 86400.0) as i64;
        let max_seconds = ((max - 2_440_587.5) * 86400.0) as i64;
        self.timestamp_subsec_nanos() < 1_000_000_000
            && self.timestamp() >= min_seconds
            && (self.timestamp() < max_seconds
                || (self.timestamp() == max_seconds && self.timestamp_subsec_nanos() == 0))
    }
}

#[cfg(feature = "alloc")]
fn start_of_date<Tz: TimeZone>(date: NaiveDate, zone: &Tz) -> Result<DateTime<Tz>> {
    let midnight = date
        .and_hms_opt(0, 0, 0)
        .ok_or_else(|| invalid_time("calendar date out of range"))?;
    // Chrono has no at-start-of-day operation. Choose the earliest occurrence,
    // advancing through a midnight gap (including an entirely skipped date).
    // Time-zone transitions have whole-second resolution.
    for seconds in 0..=86400 {
        let local = midnight
            .checked_add_signed(Duration::seconds(seconds))
            .ok_or_else(|| invalid_time("calendar date out of range"))?;
        if let Some(time) = zone.from_local_datetime(&local).earliest() {
            return Ok(time);
        }
    }
    Err(invalid_time(
        "no valid local time within a day of the date boundary",
    ))
}
