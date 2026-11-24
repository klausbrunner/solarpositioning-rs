//! Julian dates and time scales for solar positioning.
//!
//! For ΔT estimates, see [`crate::delta_t`].

#![allow(clippy::unreadable_literal)]
#![allow(clippy::many_single_char_names)]

use crate::math::floor;
use crate::{Error, Result};
#[cfg(feature = "jiff")]
use jiff::{SignedDuration, Timestamp};

/// Seconds per day (86,400)
const SECONDS_PER_DAY: f64 = 86_400.0;

#[cfg(feature = "jiff")]
const UNIX_EPOCH_JD: f64 = 2_440_587.5;

/// Julian Day Number for J2000.0 epoch (2000-01-01 12:00:00 UTC)
const J2000_JDN: f64 = 2_451_545.0;

/// Days per Julian century
const DAYS_PER_CENTURY: f64 = 36_525.0;

/// Inclusive UT and TT limits shared by positions and events.
#[derive(Debug, Clone, Copy)]
pub(crate) struct TimeRange {
    start: f64,
    end: f64,
}

impl TimeRange {
    pub(crate) const SPA: Self = Self::new(-2000, 6000);
    pub(crate) const GRENA3: Self = Self::new(2010, 2110);

    pub(crate) const fn new(first_year: i32, last_year: i32) -> Self {
        Self {
            start: year_start(first_year),
            end: year_start(last_year + 1),
        }
    }

    pub(crate) fn contains_julian(self, time: f64, delta_t: f64) -> bool {
        (self.start..=self.end).contains(&time)
            && (self.start..=self.end).contains(&(time + delta_t / SECONDS_PER_DAY))
    }

    #[cfg(feature = "jiff")]
    pub(crate) fn contains_timestamp(self, time: Timestamp, delta_t: f64) -> bool {
        let start = ((self.start - UNIX_EPOCH_JD) * SECONDS_PER_DAY) as i64;
        let end = ((self.end - UNIX_EPOCH_JD) * SECONDS_PER_DAY) as i64;
        let since_epoch = time.as_duration();
        let elapsed = since_epoch - SignedDuration::from_secs(start);
        let remaining = SignedDuration::from_secs(end) - since_epoch;
        // Compare durations to delta T before converting to Julian days, preserving boundaries.
        !elapsed.is_negative()
            && !remaining.is_negative()
            && delta_t >= -elapsed.as_secs_f64()
            && delta_t <= remaining.as_secs_f64()
    }
}

// Proleptic Gregorian January 1, including astronomical year zero.
const fn year_start(year: i32) -> f64 {
    let y = year - 1;
    1_721_425.5 + (365 * y + y.div_euclid(4) - y.div_euclid(100) + y.div_euclid(400)) as f64
}

/// Validates UTC date/time components against both field ranges and the calendar.
fn validate_utc_components(
    year: i32,
    month: u32,
    day: u32,
    hour: u32,
    minute: u32,
    second: f64,
) -> Result<()> {
    if !(1..=12).contains(&month) {
        return Err(Error::InvalidDateTime {
            message: "month must be between 1 and 12",
        });
    }
    if !(1..=31).contains(&day) {
        return Err(Error::InvalidDateTime {
            message: "day must be between 1 and 31",
        });
    }
    if hour > 23 {
        return Err(Error::InvalidDateTime {
            message: "hour must be between 0 and 23",
        });
    }
    if minute > 59 {
        return Err(Error::InvalidDateTime {
            message: "minute must be between 0 and 59",
        });
    }
    if !(0.0..60.0).contains(&second) {
        return Err(Error::InvalidDateTime {
            message: "second must be between 0 and 59.999...",
        });
    }
    if day > days_in_month(year, month) {
        return Err(Error::InvalidDateTime {
            message: "day is out of range for month",
        });
    }

    Ok(())
}

/// Julian date representation for astronomical calculations.
///
/// Supports both Julian Date (JD) and Julian Ephemeris Date (JDE) calculations.
/// Calendar constructors use the proleptic Gregorian calendar, including before 1582,
/// with year 0 representing 1 BCE. UTC approximates UT1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JulianDate {
    /// Julian Date (JD) - referenced to UT1
    jd: f64,
    /// Delta T in seconds - difference between TT and UT1
    delta_t: f64,
}

impl JulianDate {
    /// Creates a Julian date directly from continuous UT and delta T (TT minus UT1).
    ///
    /// This bypasses calendar conversion.
    ///
    /// # Errors
    /// Returns an error if either value, or the resulting TT Julian date, is non-finite.
    pub fn new(julian_date: f64, delta_t: f64) -> Result<Self> {
        if !julian_date.is_finite()
            || !delta_t.is_finite()
            || !(julian_date + delta_t / SECONDS_PER_DAY).is_finite()
        {
            return Err(Error::InvalidDateTime {
                message: "Julian date and delta_t must be finite",
            });
        }
        Ok(Self {
            jd: julian_date,
            delta_t,
        })
    }

    /// Creates a Julian date from a Jiff timestamp, approximating UT1 with UTC.
    ///
    /// `delta_t` is TT minus UT1 in seconds.
    ///
    /// # Errors
    /// Returns an error for non-finite `delta_t`.
    #[cfg(feature = "jiff")]
    #[cfg_attr(docsrs, doc(cfg(feature = "jiff")))]
    pub fn from_timestamp(time: Timestamp, delta_t: f64) -> Result<Self> {
        Self::new(timestamp_to_julian(time), delta_t)
    }

    /// Creates a new Julian date from year, month, day, hour, minute, and second in UTC.
    ///
    /// # Arguments
    /// * `year` - Year (can be negative for BCE years)
    /// * `month` - Month (1-12)
    /// * `day` - Day of month (1-31)
    /// * `hour` - Hour (0-23)
    /// * `minute` - Minute (0-59)
    /// * `second` - Second (0-59, can include fractional seconds)
    /// * `delta_t` - ΔT in seconds (difference between TT and UT1)
    ///
    /// # Returns
    /// Julian date or error if the date is invalid
    ///
    /// # Errors
    /// Returns error if any date/time component is outside valid ranges (month 1-12, day 1-31,
    /// hour 0-23, minute 0-59, second 0-59.999) or if `delta_t` is not finite.
    ///
    /// # Example
    /// ```
    /// # use solar_positioning::time::JulianDate;
    /// let jd = JulianDate::from_utc(2023, 6, 21, 12, 0, 0.0, 69.0).unwrap();
    /// assert!(jd.julian_date() > 2_460_000.0);
    /// ```
    pub fn from_utc(
        year: i32,
        month: u32,
        day: u32,
        hour: u32,
        minute: u32,
        second: f64,
        delta_t: f64,
    ) -> Result<Self> {
        validate_utc_components(year, month, day, hour, minute, second)?;
        Self::new(
            calculate_julian_date(year, month, day, hour, minute, second),
            delta_t,
        )
    }

    /// Creates a Julian date assuming ΔT = 0.
    ///
    /// # Arguments
    /// * `year` - Year (can be negative for BCE years)
    /// * `month` - Month (1-12)
    /// * `day` - Day of month (1-31)
    /// * `hour` - Hour (0-23)
    /// * `minute` - Minute (0-59)
    /// * `second` - Second (0-59, can include fractional seconds)
    ///
    /// # Returns
    /// Returns `Ok(JulianDate)` with ΔT = 0 on success.
    ///
    /// # Errors
    /// Returns error if the date/time components are outside valid ranges.
    pub fn from_utc_simple(
        year: i32,
        month: u32,
        day: u32,
        hour: u32,
        minute: u32,
        second: f64,
    ) -> Result<Self> {
        Self::from_utc(year, month, day, hour, minute, second, 0.0)
    }

    /// Gets the Julian Date (JD) value.
    ///
    /// # Returns
    /// Julian Date referenced to UT1
    #[must_use]
    pub const fn julian_date(&self) -> f64 {
        self.jd
    }

    /// Gets the ΔT value in seconds.
    ///
    /// # Returns
    /// ΔT (Delta T) in seconds
    #[must_use]
    pub const fn delta_t(&self) -> f64 {
        self.delta_t
    }

    /// Calculates the Julian Ephemeris Day (JDE).
    ///
    /// JDE = JD + ΔT/86400
    ///
    /// # Returns
    /// Julian Ephemeris Day
    #[must_use]
    pub fn julian_ephemeris_day(&self) -> f64 {
        self.jd + self.delta_t / SECONDS_PER_DAY
    }

    /// Calculates the Julian Century (JC) from J2000.0.
    ///
    /// JC = (JD - 2451545.0) / 36525
    ///
    /// # Returns
    /// Julian centuries since J2000.0 epoch
    #[must_use]
    pub fn julian_century(&self) -> f64 {
        (self.jd - J2000_JDN) / DAYS_PER_CENTURY
    }

    /// Calculates the Julian Ephemeris Century (JCE) from J2000.0.
    ///
    /// JCE = (JDE - 2451545.0) / 36525
    ///
    /// # Returns
    /// Julian ephemeris centuries since J2000.0 epoch
    #[must_use]
    pub fn julian_ephemeris_century(&self) -> f64 {
        (self.julian_ephemeris_day() - J2000_JDN) / DAYS_PER_CENTURY
    }

    /// Calculates the Julian Ephemeris Millennium (JME) from J2000.0.
    ///
    /// JME = JCE / 10
    ///
    /// # Returns
    /// Julian ephemeris millennia since J2000.0 epoch
    #[must_use]
    pub fn julian_ephemeris_millennium(&self) -> f64 {
        self.julian_ephemeris_century() / 10.0
    }
}

#[cfg(feature = "jiff")]
pub(crate) fn timestamp_to_julian(time: Timestamp) -> f64 {
    // Jiff returns signed subseconds, so floor the second for pre-epoch times.
    let nanos = time.subsec_nanosecond();
    let seconds = time.as_second() - i64::from(nanos < 0);
    // Add the fractional day as a unit so rounding cannot reverse adjacent instants.
    let midnight = UNIX_EPOCH_JD + seconds.div_euclid(86_400) as f64;
    let second_of_day =
        seconds.rem_euclid(86_400) as f64 + f64::from(nanos.rem_euclid(1_000_000_000)) / 1e9;
    midnight + second_of_day / SECONDS_PER_DAY
}

/// Meeus's Julian date calculation, applying Gregorian rules to every year.
fn calculate_julian_date(
    year: i32,
    month: u32,
    day: u32,
    hour: u32,
    minute: u32,
    second: f64,
) -> f64 {
    let mut y = f64::from(year);

    // Adjust for January and February being treated as months 13 and 14 of previous year
    let m = if month < 3 {
        y -= 1.0;
        month + 12
    } else {
        month
    };

    // Calculate fractional day
    let d = f64::from(day) + (f64::from(hour) + (f64::from(minute) + second / 60.0) / 60.0) / 24.0;

    // Basic Julian Date calculation
    let jd = floor(365.25 * (y + 4716.0)) + floor(30.6001 * f64::from(m + 1)) + d - 1524.5;
    let a = floor(y / 100.0);
    let b = 2.0 - a + floor(a / 4.0);
    jd + b
}

fn days_in_month(year: i32, month: u32) -> u32 {
    let is_leap_year = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;

    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if is_leap_year => 29,
        2 => 28,
        _ => unreachable!("month already validated"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_julian_date_invalid_day_validation() {
        assert!(JulianDate::from_utc(2024, 2, 30, 0, 0, 0.0, 0.0).is_err());
        assert!(JulianDate::from_utc(2024, 2, 29, 0, 0, 0.0, 0.0).is_ok());
        assert!(JulianDate::from_utc(1900, 2, 29, 0, 0, 0.0, 0.0).is_err());
        assert!(JulianDate::from_utc(1500, 2, 29, 0, 0, 0.0, 0.0).is_err());
        assert!(JulianDate::from_utc(0, 2, 29, 0, 0, 0.0, 0.0).is_ok());
        assert!(JulianDate::from_utc(-100, 2, 29, 0, 0, 0.0, 0.0).is_err());
    }

    #[test]
    fn test_julian_date_validation() {
        assert!(JulianDate::from_utc(2024, 13, 1, 0, 0, 0.0, 0.0).is_err()); // Invalid month
        assert!(JulianDate::from_utc(2024, 1, 32, 0, 0, 0.0, 0.0).is_err()); // Invalid day
        assert!(JulianDate::from_utc(2024, 1, 1, 24, 0, 0.0, 0.0).is_err()); // Invalid hour
        assert!(JulianDate::from_utc(2024, 1, 1, 0, 60, 0.0, 0.0).is_err()); // Invalid minute
        assert!(JulianDate::from_utc(2024, 1, 1, 0, 0, 60.0, 0.0).is_err()); // Invalid second
        assert!(JulianDate::from_utc(2024, 1, 1, 0, 0, 0.0, f64::NAN).is_err()); // Non-finite delta_t
        assert!(JulianDate::from_utc(2024, 1, 1, 0, 0, 0.0, f64::INFINITY).is_err());
        // Date conversion must not overflow at the i32 year boundaries.
        for year in [i32::MIN, i32::MAX] {
            assert!(
                JulianDate::from_utc(year, 1, 1, 0, 0, 0.0, 0.0)
                    .unwrap()
                    .julian_date()
                    .is_finite()
            );
        }
    }

    #[test]
    fn test_julian_centuries() {
        let jd = JulianDate::from_utc(2000, 1, 1, 12, 0, 0.0, 0.0).unwrap();

        // J2000.0 should give JC = 0
        assert!(jd.julian_century().abs() < EPSILON);
        assert!(jd.julian_ephemeris_century().abs() < EPSILON);
        assert!(jd.julian_ephemeris_millennium().abs() < EPSILON);
    }

    #[test]
    fn test_julian_ephemeris_day() {
        let delta_t = 69.0; // seconds
        let jd = JulianDate::from_utc(2023, 6, 21, 12, 0, 0.0, delta_t).unwrap();

        let jde = jd.julian_ephemeris_day();
        let expected = jd.julian_date() + delta_t / SECONDS_PER_DAY;

        assert!((jde - expected).abs() < EPSILON);
    }

    #[test]
    fn test_proleptic_gregorian_dates() {
        // Shared with Java; Julian days cross-check against JDK JulianFields.
        for (year, month, day, hour, expected) in [
            (-4713, 11, 24, 0, -0.5),
            (-123, 12, 28, 0, 1_676_496.5),
            (-123, 12, 29, 0, 1_676_497.5),
            (837, 4, 14, 0, 2_026_871.5),
            (1582, 10, 4, 0, 2_299_149.5),
            (1582, 10, 10, 0, 2_299_155.5),
            (1582, 10, 15, 0, 2_299_160.5),
            (1970, 1, 1, 0, 2_440_587.5),
            (2000, 1, 1, 12, J2000_JDN),
        ] {
            let jd = JulianDate::from_utc(year, month, day, hour, 0, 0.0, 69.184).unwrap();
            assert_eq!(jd.julian_date(), expected, "{year}-{month}-{day}");
            assert_eq!(jd.delta_t(), 69.184);
            #[cfg(feature = "jiff")]
            {
                let time = jiff::tz::Offset::UTC
                    .to_timestamp(jiff::civil::datetime(
                        year as i16,
                        month as i8,
                        day as i8,
                        hour as i8,
                        0,
                        0,
                        0,
                    ))
                    .unwrap();
                assert_eq!(JulianDate::from_timestamp(time, 69.184).unwrap(), jd);
            }
        }
    }

    #[cfg(feature = "jiff")]
    #[test]
    fn timestamp_to_julian_is_monotone_across_second_boundaries() {
        for seconds in [-90_000_000_024, 0, 90_000_000_024] {
            let after = Timestamp::from_second(seconds).unwrap();
            let before = after - SignedDuration::from_nanos(1);
            assert!(timestamp_to_julian(before) <= timestamp_to_julian(after));
        }
    }

    #[cfg(feature = "jiff")]
    #[test]
    fn timestamp_matches_numeric_fractional_seconds() {
        for year in [-2000, 0, 1500, 1582, 1969, 1970, 2024, 6000] {
            let time = jiff::tz::Offset::UTC
                .to_timestamp(jiff::civil::datetime(year, 12, 31, 23, 59, 59, 123_456_789))
                .unwrap();
            let numeric =
                JulianDate::from_utc(i32::from(year), 12, 31, 23, 59, 59.123_456_789, 69.184)
                    .unwrap();
            let from_timestamp = JulianDate::from_timestamp(time, 69.184).unwrap();
            // The different arithmetic paths can round by one Julian-day ulp.
            assert!((numeric.julian_date() - from_timestamp.julian_date()).abs() < 1e-9);
        }
    }
}
