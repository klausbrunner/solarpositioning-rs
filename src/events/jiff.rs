#[cfg(feature = "alloc")]
use super::invalid_time;
use super::{EventPosition, SearchTime, SolarEvents};
use crate::{
    Horizon, Location, Result,
    time::{JulianDate, TimeRange, timestamp_to_julian},
};
use ::jiff::{SignedDuration, Timestamp};
#[cfg(feature = "alloc")]
use {
    super::Events,
    ::jiff::{
        ZonedDateTime,
        civil::Date,
        tz::{AmbiguousOffset, TimeZone},
    },
    alloc::vec::Vec,
    core::ops::Range,
};

impl<P> SolarEvents<P>
where
    P: Fn(JulianDate, Location) -> Result<EventPosition>,
{
    /// Finds the next rising crossing in `(start, end]`.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_rise_from_julian`].
    pub fn next_rise(
        &self,
        start: Timestamp,
        end: Timestamp,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Option<Timestamp>> {
        self.crossing(start, end, location, delta_t, horizon, 1.0)
    }

    /// Finds the next setting crossing in `(start, end]`.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_rise_from_julian`].
    pub fn next_set(
        &self,
        start: Timestamp,
        end: Timestamp,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Option<Timestamp>> {
        self.crossing(start, end, location, delta_t, horizon, -1.0)
    }

    /// Finds the next upper meridian transit in `(start, end]`.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_transit_from_julian`].
    pub fn next_transit(
        &self,
        start: Timestamp,
        end: Timestamp,
        longitude: f64,
        delta_t: f64,
    ) -> Result<Option<Timestamp>> {
        self.transit(start, end, longitude, delta_t)
    }

    /// Collects all events in a local calendar date, including its start and excluding its end.
    ///
    /// Handles daylight-saving changes, repeated dates and skipped dates. Returned
    /// times use `zone`. Pass [`Horizon::SunriseSunset`] for standard rise/set times.
    /// Historical rollbacks that split a local date into separate intervals are not supported.
    ///
    /// # Errors
    /// Has the same validation and provider errors as [`Self::next_rise_from_julian`].
    /// Also returns an error if a calendar date boundary cannot be resolved.
    #[cfg(feature = "alloc")]
    #[cfg_attr(docsrs, doc(cfg(all(feature = "jiff", feature = "alloc"))))]
    pub fn for_date(
        &self,
        date: Date,
        zone: &TimeZone,
        location: Location,
        delta_t: f64,
        horizon: Horizon,
    ) -> Result<Events<ZonedDateTime>> {
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
    #[cfg_attr(docsrs, doc(cfg(all(feature = "jiff", feature = "alloc"))))]
    pub fn for_date_multiple(
        &self,
        date: Date,
        zone: &TimeZone,
        location: Location,
        delta_t: f64,
        horizons: impl IntoIterator<Item = Horizon>,
    ) -> Result<Vec<(Horizon, Events<ZonedDateTime>)>> {
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
impl Events<Timestamp> {
    fn into_timezone(self, zone: &TimeZone) -> Events<ZonedDateTime> {
        let convert = |time: Timestamp| time.to_zoned_datetime(zone.clone());
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
fn date_interval(date: Date, zone: &TimeZone) -> Result<Range<Timestamp>> {
    let next_date = date
        .tomorrow()
        .map_err(|_| invalid_time("calendar date out of range"))?;
    Ok(start_of_date(date, zone)?..start_of_date(next_date, zone)?)
}

impl SearchTime for Timestamp {
    fn julian_day(self) -> f64 {
        timestamp_to_julian(self)
    }

    fn hours_until(self, end: Self) -> f64 {
        end.duration_since(self).as_secs_f64() / 3600.0
    }

    fn shift_hours(self, hours: f64) -> Self {
        self + SignedDuration::from_secs_f64(hours * 3600.0)
    }

    fn in_range(self, range: TimeRange, delta_t: f64) -> bool {
        range.contains_timestamp(self, delta_t)
    }
}

#[cfg(feature = "alloc")]
fn start_of_date(date: Date, zone: &TimeZone) -> Result<Timestamp> {
    let midnight = zone.to_ambiguous_timestamp(date.at(0, 0, 0, 0));
    let earlier = midnight
        .earlier()
        .map_err(|_| invalid_time("calendar date out of range"))?;
    if matches!(midnight.offset(), AmbiguousOffset::Gap { .. }) {
        // Use the transition itself: shifting midnight by the gap length can
        // skip valid time when the gap starts before midnight.
        return zone
            .following(earlier)
            .next()
            .map(|transition| transition.timestamp())
            .ok_or_else(|| invalid_time("no time-zone transition at date boundary"));
    }
    Ok(earlier)
}
