//! Time-related calculations for solar positioning.
//!
//! This module provides Julian date calculations and ΔT (Delta T) estimation
//! following the algorithms from NREL SPA and Espenak & Meeus.

#![allow(clippy::unreadable_literal)]
#![allow(clippy::many_single_char_names)]

use crate::math::{floor, polynomial};
use crate::{Error, Result};
#[cfg(feature = "chrono")]
use chrono::{Datelike, TimeZone, Timelike};

/// Seconds per day (86,400)
const SECONDS_PER_DAY: f64 = 86_400.0;

/// Julian Day Number for J2000.0 epoch (2000-01-01 12:00:00 UTC)
const J2000_JDN: f64 = 2_451_545.0;

/// Days per Julian century
const DAYS_PER_CENTURY: f64 = 36_525.0;

/// Julian date representation for astronomical calculations.
///
/// Follows the SPA algorithm described in Reda & Andreas (2003).
/// Supports both Julian Date (JD) and Julian Ephemeris Date (JDE) calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JulianDate {
    /// Julian Date (JD) - referenced to UT1
    jd: f64,
    /// Delta T in seconds - difference between TT and UT1
    delta_t: f64,
}

impl JulianDate {
    /// Creates a new Julian date from a timezone-aware chrono `DateTime`.
    ///
    /// Converts datetime to UTC for proper Julian Date calculation.
    ///
    /// # Arguments
    /// * `datetime` - Timezone-aware date and time
    /// * `delta_t` - ΔT in seconds (difference between TT and UT1)
    ///
    /// # Returns
    /// Returns `Ok(JulianDate)` on success.
    ///
    /// # Errors
    /// Returns error if the date/time components are invalid (e.g., invalid month, day, hour).
    #[cfg(feature = "chrono")]
    pub fn from_datetime<Tz: TimeZone>(
        datetime: &chrono::DateTime<Tz>,
        delta_t: f64,
    ) -> Result<Self> {
        // Convert the entire datetime to UTC for proper Julian Date calculation
        let utc_datetime = datetime.with_timezone(&chrono::Utc);
        Self::from_utc(
            utc_datetime.year(),
            utc_datetime.month(),
            utc_datetime.day(),
            utc_datetime.hour(),
            utc_datetime.minute(),
            f64::from(utc_datetime.second()) + f64::from(utc_datetime.nanosecond()) / 1e9,
            delta_t,
        )
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
    /// Returns error if any date/time component is outside valid ranges (month 1-12, day 1-31, hour 0-23, minute 0-59, second 0-59.999).
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
        // Validate input ranges
        if !(1..=12).contains(&month) {
            return Err(Error::invalid_datetime("month must be between 1 and 12"));
        }
        if !(1..=31).contains(&day) {
            return Err(Error::invalid_datetime("day must be between 1 and 31"));
        }
        if hour > 23 {
            return Err(Error::invalid_datetime("hour must be between 0 and 23"));
        }
        if minute > 59 {
            return Err(Error::invalid_datetime("minute must be between 0 and 59"));
        }
        if !(0.0..60.0).contains(&second) {
            return Err(Error::invalid_datetime(
                "second must be between 0 and 59.999...",
            ));
        }

        if day > days_in_month(year, month, day)? {
            return Err(Error::invalid_datetime("day is out of range for month"));
        }

        let jd = calculate_julian_date(year, month, day, hour, minute, second);
        Ok(Self { jd, delta_t })
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

    /// Add days to the Julian date (like Java constructor: new `JulianDate(jd.julianDate()` + i - 1, 0))
    pub(crate) fn add_days(self, days: f64) -> Self {
        Self {
            jd: self.jd + days,
            delta_t: self.delta_t,
        }
    }
}

/// Calculates Julian Date from UTC date/time components.
///
/// This follows the algorithm from Reda & Andreas (2003), which is based on
/// Meeus, "Astronomical Algorithms", 2nd edition.
fn calculate_julian_date(
    year: i32,
    month: u32,
    day: u32,
    hour: u32,
    minute: u32,
    second: f64,
) -> f64 {
    let mut y = year;
    let mut m = i32::try_from(month).expect("month should be valid i32");

    // Adjust for January and February being treated as months 13 and 14 of previous year
    if m < 3 {
        y -= 1;
        m += 12;
    }

    // Calculate fractional day
    let d = f64::from(day) + (f64::from(hour) + (f64::from(minute) + second / 60.0) / 60.0) / 24.0;

    // Basic Julian Date calculation
    let mut jd =
        floor(365.25 * (f64::from(y) + 4716.0)) + floor(30.6001 * f64::from(m + 1)) + d - 1524.5;

    // Gregorian calendar correction (after October 15, 1582)
    // JDN 2299161 corresponds to October 15, 1582
    if jd >= 2_299_161.0 {
        let a = floor(f64::from(y) / 100.0);
        let b = 2.0 - a + floor(a / 4.0);
        jd += b;
    }

    jd
}

const fn is_gregorian_date(year: i32, month: u32, day: u32) -> bool {
    year > 1582 || (year == 1582 && (month > 10 || (month == 10 && day >= 15)))
}

const fn is_leap_year(year: i32, is_gregorian: bool) -> bool {
    if is_gregorian {
        (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
    } else {
        year % 4 == 0
    }
}

fn days_in_month(year: i32, month: u32, day: u32) -> Result<u32> {
    if year == 1582 && month == 10 && (5..=14).contains(&day) {
        return Err(Error::invalid_datetime(
            "dates 1582-10-05 through 1582-10-14 do not exist in Gregorian calendar",
        ));
    }

    let is_gregorian = is_gregorian_date(year, month, day);
    let days = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if is_leap_year(year, is_gregorian) {
                29
            } else {
                28
            }
        }
        _ => unreachable!("month already validated"),
    };
    Ok(days)
}

/// ΔT (Delta T) estimation functions.
///
/// ΔT represents the difference between Terrestrial Time (TT) and Universal Time (UT1).
/// These estimates are based on Espenak and Meeus polynomial fits updated in 2014.
pub struct DeltaT;

impl DeltaT {
    /// Estimates ΔT for a given decimal year.
    ///
    /// Based on polynomial fits from Espenak & Meeus, updated 2014.
    /// See: <https://www.eclipsewise.com/help/deltatpoly2014.html>
    ///
    /// # Arguments
    /// * `decimal_year` - Year with fractional part (e.g., 2024.5 for mid-2024)
    ///
    /// # Returns
    /// Estimated ΔT in seconds
    ///
    /// # Errors
    /// Returns error for years outside the valid range (-500 to 3000 CE)
    ///
    /// # Example
    /// ```
    /// # use solar_positioning::time::DeltaT;
    /// let delta_t = DeltaT::estimate(2024.0).unwrap();
    /// assert!(delta_t > 60.0 && delta_t < 80.0); // Reasonable range for 2024
    /// ```
    #[allow(clippy::too_many_lines)] // Comprehensive polynomial fit across historical periods
    pub fn estimate(decimal_year: f64) -> Result<f64> {
        let year = decimal_year;

        if !year.is_finite() {
            return Err(Error::invalid_datetime("year must be finite"));
        }

        let delta_t = if year < -500.0 {
            let u = (year - 1820.0) / 100.0;
            polynomial(&[-20.0, 0.0, 32.0], u)
        } else if year < 500.0 {
            let u = year / 100.0;
            polynomial(
                &[
                    10583.6,
                    -1014.41,
                    33.78311,
                    -5.952053,
                    -0.1798452,
                    0.022174192,
                    0.0090316521,
                ],
                u,
            )
        } else if year < 1600.0 {
            let u = (year - 1000.0) / 100.0;
            polynomial(
                &[
                    1574.2,
                    -556.01,
                    71.23472,
                    0.319781,
                    -0.8503463,
                    -0.005050998,
                    0.0083572073,
                ],
                u,
            )
        } else if year < 1700.0 {
            let t = year - 1600.0;
            polynomial(&[120.0, -0.9808, -0.01532, 1.0 / 7129.0], t)
        } else if year < 1800.0 {
            let t = year - 1700.0;
            polynomial(
                &[8.83, 0.1603, -0.0059285, 0.00013336, -1.0 / 1_174_000.0],
                t,
            )
        } else if year < 1860.0 {
            let t = year - 1800.0;
            polynomial(
                &[
                    13.72,
                    -0.332447,
                    0.0068612,
                    0.0041116,
                    -0.00037436,
                    0.0000121272,
                    -0.0000001699,
                    0.000000000875,
                ],
                t,
            )
        } else if year < 1900.0 {
            let t = year - 1860.0;
            polynomial(
                &[
                    7.62,
                    0.5737,
                    -0.251754,
                    0.01680668,
                    -0.0004473624,
                    1.0 / 233_174.0,
                ],
                t,
            )
        } else if year < 1920.0 {
            let t = year - 1900.0;
            polynomial(&[-2.79, 1.494119, -0.0598939, 0.0061966, -0.000197], t)
        } else if year < 1941.0 {
            let t = year - 1920.0;
            polynomial(&[21.20, 0.84493, -0.076100, 0.0020936], t)
        } else if year < 1961.0 {
            let t = year - 1950.0;
            polynomial(&[29.07, 0.407, -1.0 / 233.0, 1.0 / 2547.0], t)
        } else if year < 1986.0 {
            let t = year - 1975.0;
            polynomial(&[45.45, 1.067, -1.0 / 260.0, -1.0 / 718.0], t)
        } else if year < 2005.0 {
            let t = year - 2000.0;
            polynomial(
                &[
                    63.86,
                    0.3345,
                    -0.060374,
                    0.0017275,
                    0.000651814,
                    0.00002373599,
                ],
                t,
            )
        } else if year < 2015.0 {
            let t = year - 2005.0;
            polynomial(&[64.69, 0.2930], t)
        } else if year <= 3000.0 {
            let t = year - 2015.0;
            polynomial(&[67.62, 0.3645, 0.0039755], t)
        } else {
            return Err(Error::invalid_datetime(
                "ΔT estimates not available beyond year 3000",
            ));
        };

        Ok(delta_t)
    }

    /// Estimates ΔT from year and month.
    ///
    /// Calculates decimal year as: year + (month - 0.5) / 12
    ///
    /// # Arguments
    /// * `year` - Year
    /// * `month` - Month (1-12)
    ///
    /// # Returns
    /// Returns estimated ΔT in seconds.
    ///
    /// # Errors
    /// Returns error if month is outside the range 1-12.
    ///
    /// # Panics
    /// This function does not panic.
    pub fn estimate_from_date(year: i32, month: u32) -> Result<f64> {
        if !(1..=12).contains(&month) {
            return Err(Error::invalid_datetime("month must be between 1 and 12"));
        }

        let decimal_year = f64::from(year) + (f64::from(month) - 0.5) / 12.0;
        Self::estimate(decimal_year)
    }

    /// Estimates ΔT from any date-like type.
    ///
    /// Convenience method that extracts the year and month from any chrono type
    /// that implements `Datelike` (`DateTime`, `NaiveDateTime`, `NaiveDate`, etc.).
    ///
    /// # Arguments
    /// * `date` - Any date-like type
    ///
    /// # Returns
    /// Returns estimated ΔT in seconds.
    ///
    /// # Errors
    /// Returns error if the date components are invalid.
    ///
    /// # Example
    /// ```
    /// # use solar_positioning::time::DeltaT;
    /// # use chrono::{DateTime, FixedOffset, NaiveDate};
    ///
    /// // Works with DateTime
    /// let datetime = "2024-06-21T12:00:00-07:00".parse::<DateTime<FixedOffset>>().unwrap();
    /// let delta_t = DeltaT::estimate_from_date_like(datetime).unwrap();
    /// assert!(delta_t > 60.0 && delta_t < 80.0);
    ///
    /// // Also works with NaiveDate
    /// let date = NaiveDate::from_ymd_opt(2024, 6, 21).unwrap();
    /// let delta_t2 = DeltaT::estimate_from_date_like(date).unwrap();
    /// assert_eq!(delta_t, delta_t2);
    #[cfg(feature = "chrono")]
    #[allow(clippy::needless_pass_by_value)]
    pub fn estimate_from_date_like<D: Datelike>(date: D) -> Result<f64> {
        Self::estimate_from_date(date.year(), date.month())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_julian_date_creation() {
        let jd = JulianDate::from_utc(2000, 1, 1, 12, 0, 0.0, 0.0).unwrap();

        // J2000.0 epoch should be exactly 2451545.0
        assert!((jd.julian_date() - J2000_JDN).abs() < EPSILON);
        assert_eq!(jd.delta_t(), 0.0);
    }

    #[test]
    fn test_julian_date_invalid_day_validation() {
        assert!(JulianDate::from_utc(2024, 2, 30, 0, 0, 0.0, 0.0).is_err());
        assert!(JulianDate::from_utc(2024, 2, 29, 0, 0, 0.0, 0.0).is_ok());
        assert!(JulianDate::from_utc(1900, 2, 29, 0, 0, 0.0, 0.0).is_err());
        assert!(JulianDate::from_utc(1500, 2, 29, 0, 0, 0.0, 0.0).is_ok());
        assert!(JulianDate::from_utc(1582, 10, 10, 0, 0, 0.0, 0.0).is_err());
        assert!(JulianDate::from_utc(1582, 10, 4, 0, 0, 0.0, 0.0).is_ok());
        assert!(JulianDate::from_utc(1582, 10, 15, 0, 0, 0.0, 0.0).is_ok());
    }

    #[test]
    fn test_julian_date_validation() {
        assert!(JulianDate::from_utc(2024, 13, 1, 0, 0, 0.0, 0.0).is_err()); // Invalid month
        assert!(JulianDate::from_utc(2024, 1, 32, 0, 0, 0.0, 0.0).is_err()); // Invalid day
        assert!(JulianDate::from_utc(2024, 1, 1, 24, 0, 0.0, 0.0).is_err()); // Invalid hour
        assert!(JulianDate::from_utc(2024, 1, 1, 0, 60, 0.0, 0.0).is_err()); // Invalid minute
        assert!(JulianDate::from_utc(2024, 1, 1, 0, 0, 60.0, 0.0).is_err()); // Invalid second
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
    fn test_gregorian_calendar_correction() {
        // Test dates before and after Gregorian calendar adoption
        // October 4, 1582 was followed by October 15, 1582
        let julian_date = JulianDate::from_utc(1582, 10, 4, 12, 0, 0.0, 0.0).unwrap();
        let gregorian_date = JulianDate::from_utc(1582, 10, 15, 12, 0, 0.0, 0.0).unwrap();

        // The calendar dates are 11 days apart, but in Julian Day Numbers they should be 1 day apart
        // because the 10-day gap was artificial
        let diff = gregorian_date.julian_date() - julian_date.julian_date();
        assert!(
            (diff - 1.0).abs() < 1e-6,
            "Expected 1 day difference in JD, got {diff}"
        );

        // Test that the Gregorian correction is applied correctly
        // Dates after October 15, 1582 should have the correction
        let pre_gregorian = JulianDate::from_utc(1582, 10, 1, 12, 0, 0.0, 0.0).unwrap();
        let post_gregorian = JulianDate::from_utc(1583, 1, 1, 12, 0, 0.0, 0.0).unwrap();

        // Verify that both exist and the calculation doesn't panic
        assert!(pre_gregorian.julian_date() > 2_000_000.0);
        assert!(post_gregorian.julian_date() > pre_gregorian.julian_date());
    }

    #[test]
    fn test_delta_t_modern_estimates() {
        // Test some known ranges
        let delta_t_2000 = DeltaT::estimate(2000.0).unwrap();
        let delta_t_2020 = DeltaT::estimate(2020.0).unwrap();

        assert!(delta_t_2000 > 60.0 && delta_t_2000 < 70.0);
        assert!(delta_t_2020 > 65.0 && delta_t_2020 < 75.0);
        assert!(delta_t_2020 > delta_t_2000); // ΔT is generally increasing
    }

    #[test]
    fn test_delta_t_historical_estimates() {
        let delta_t_1900 = DeltaT::estimate(1900.0).unwrap();
        let delta_t_1950 = DeltaT::estimate(1950.0).unwrap();

        assert!(delta_t_1900 < 0.0); // Negative in early 20th century
        assert!(delta_t_1950 > 25.0 && delta_t_1950 < 35.0);
    }

    #[test]
    fn test_delta_t_boundary_conditions() {
        // Test edge cases
        assert!(DeltaT::estimate(-500.0).is_ok());
        assert!(DeltaT::estimate(3000.0).is_ok());
        assert!(DeltaT::estimate(-501.0).is_ok()); // Should work for ancient dates
        assert!(DeltaT::estimate(3001.0).is_err()); // Should fail beyond 3000
    }

    #[test]
    fn test_delta_t_from_date() {
        let delta_t = DeltaT::estimate_from_date(2024, 6).unwrap();
        let delta_t_decimal = DeltaT::estimate(2024.5 - 1.0 / 24.0).unwrap(); // June = month 6, so (6-0.5)/12 ≈ 0.458

        // Should be very close
        assert!((delta_t - delta_t_decimal).abs() < 0.01);

        // Test invalid month
        assert!(DeltaT::estimate_from_date(2024, 13).is_err());
        assert!(DeltaT::estimate_from_date(2024, 0).is_err());
    }

    #[test]
    fn test_delta_t_from_date_like() {
        use chrono::{DateTime, FixedOffset, NaiveDate, Utc};

        // Test with DateTime<FixedOffset>
        let datetime_fixed = "2024-06-15T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let delta_t_fixed = DeltaT::estimate_from_date_like(datetime_fixed).unwrap();

        // Test with DateTime<Utc>
        let datetime_utc = "2024-06-15T19:00:00Z".parse::<DateTime<Utc>>().unwrap();
        let delta_t_utc = DeltaT::estimate_from_date_like(datetime_utc).unwrap();

        // Test with NaiveDate
        let naive_date = NaiveDate::from_ymd_opt(2024, 6, 15).unwrap();
        let delta_t_naive_date = DeltaT::estimate_from_date_like(naive_date).unwrap();

        // Test with NaiveDateTime
        let naive_datetime = naive_date.and_hms_opt(12, 0, 0).unwrap();
        let delta_t_naive_datetime = DeltaT::estimate_from_date_like(naive_datetime).unwrap();

        // Should all be identical since we only use year/month
        assert_eq!(delta_t_fixed, delta_t_utc);
        assert_eq!(delta_t_fixed, delta_t_naive_date);
        assert_eq!(delta_t_fixed, delta_t_naive_datetime);

        // Should match estimate_from_date
        let delta_t_date = DeltaT::estimate_from_date(2024, 6).unwrap();
        assert_eq!(delta_t_fixed, delta_t_date);

        // Verify reasonable range for 2024
        assert!(delta_t_fixed > 60.0 && delta_t_fixed < 80.0);
    }

    #[test]
    fn test_specific_julian_dates() {
        // Test some well-known dates

        // Unix epoch: 1970-01-01 00:00:00 UTC
        let unix_epoch = JulianDate::from_utc(1970, 1, 1, 0, 0, 0.0, 0.0).unwrap();
        assert!((unix_epoch.julian_date() - 2_440_587.5).abs() < 1e-6);

        // Y2K: 2000-01-01 00:00:00 UTC
        let y2k = JulianDate::from_utc(2000, 1, 1, 0, 0, 0.0, 0.0).unwrap();
        assert!((y2k.julian_date() - 2_451_544.5).abs() < 1e-6);
    }
}
