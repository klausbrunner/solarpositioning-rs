//! ΔT (Delta T) estimation functions.
//!
//! ΔT represents the difference between Terrestrial Time (TT) and Universal Time (UT1).
//! Based on Espenak and Meeus, with [Espenak's 2014 update](https://www.eclipsewise.com/help/deltatpoly2014.html).
//! The branches from 2015 onwards use a [2026 adaptation](https://klaus.brunners.name/posts/delta-t-polynomials/):
//! a quartic fitted to IERS observations through mid-2026, followed by a quadratic
//! fitted to a median of simulations by Agnew (2026).
//!
//! Future values are (very) uncertain extrapolations, particularly beyond 2100.

#![allow(clippy::unreadable_literal)]

use crate::math::polynomial;
use crate::{Error, Result};
#[cfg(feature = "chrono")]
use chrono::Datelike;

/// Estimates ΔT for a given decimal year.
///
/// Uses the historical polynomials and updated fits described in this module.
///
/// Returns seconds. The year includes its fractional part (e.g. 2024.5 for mid-2024).
///
/// # Errors
/// Returns an error for non-finite years or years outside -500 through 3000.
///
/// # Example
/// ```
/// # use solar_positioning::delta_t;
/// let delta_t = delta_t::estimate(2024.0).unwrap();
/// assert!(delta_t > 60.0 && delta_t < 80.0); // Reasonable range for 2024
/// ```
#[allow(clippy::too_many_lines)] // Comprehensive polynomial fit across historical periods
pub fn estimate(decimal_year: f64) -> Result<f64> {
    let year = decimal_year;

    if !year.is_finite() {
        return Err(Error::InvalidDateTime {
            message: "year must be finite",
        });
    }

    if year < -500.0 {
        return Err(Error::InvalidDateTime {
            message: "ΔT estimates not available before year -500",
        });
    }

    let delta_t = if year < 500.0 {
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
    } else if year < 2026.5 {
        let t = year - 2015.0;
        // Retain full precision; the branches join in value and slope at 2015 and 2026.5.
        polynomial(
            &[
                67.62,
                0.2930,
                0.08753166427153103,
                -0.020049795884351372,
                0.0009758496126416308,
            ],
            t,
        )
    } else if year <= 3000.0 {
        let t = year - 2026.5;
        polynomial(&[69.14, 0.28805287963416737, 0.0057380400318460655], t)
    } else {
        return Err(Error::InvalidDateTime {
            message: "ΔT estimates not available beyond year 3000",
        });
    };

    Ok(delta_t)
}

/// Estimates ΔT at the midpoint of the given calendar month.
///
/// Calculates decimal year as: year + (month - 0.5) / 12
///
/// Returns seconds. Month is 1 through 12.
///
/// # Errors
/// Returns an error if the month is invalid or its midpoint falls outside
/// the year range supported by [`estimate`].
pub fn estimate_from_date(year: i32, month: u32) -> Result<f64> {
    if !(1..=12).contains(&month) {
        return Err(Error::InvalidDateTime {
            message: "month must be between 1 and 12",
        });
    }

    let decimal_year = f64::from(year) + (f64::from(month) - 0.5) / 12.0;
    estimate(decimal_year)
}

/// Estimates ΔT from any date-like type.
///
/// Extracts the year and month from any chrono type
/// that implements `Datelike` (`DateTime`, `NaiveDateTime`, `NaiveDate`, etc.).
/// Uses the midpoint of its calendar month, ignoring the day and time.
///
/// Returns seconds.
///
/// # Errors
/// Has the same validation errors as [`estimate_from_date`].
///
/// # Example
/// ```
/// # use solar_positioning::delta_t;
/// # use chrono::{DateTime, FixedOffset, NaiveDate};
///
/// // Works with DateTime
/// let datetime = "2024-06-21T12:00:00-07:00".parse::<DateTime<FixedOffset>>().unwrap();
/// let delta_t = delta_t::estimate_from_date_like(datetime).unwrap();
/// assert!(delta_t > 60.0 && delta_t < 80.0);
///
/// // Also works with NaiveDate
/// let date = NaiveDate::from_ymd_opt(2024, 6, 21).unwrap();
/// let delta_t2 = delta_t::estimate_from_date_like(date).unwrap();
/// assert_eq!(delta_t, delta_t2);
/// ```
#[cfg(feature = "chrono")]
#[cfg_attr(docsrs, doc(cfg(feature = "chrono")))]
#[allow(clippy::needless_pass_by_value)]
pub fn estimate_from_date_like<D: Datelike>(date: D) -> Result<f64> {
    estimate_from_date(date.year(), date.month())
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_delta_t_modern_estimates() {
        // Test some known ranges
        let delta_t_2000 = estimate(2000.0).unwrap();
        let delta_t_2020 = estimate(2020.0).unwrap();

        assert!(delta_t_2000 > 60.0 && delta_t_2000 < 70.0);
        assert!(delta_t_2020 > 65.0 && delta_t_2020 < 75.0);
        assert!(delta_t_2020 > delta_t_2000); // ΔT is generally increasing
    }

    #[test]
    fn test_delta_t_adaptation_reference_values() {
        let cases = [
            (2000.0, 63.86),
            (2005.0, 64.69),
            (2010.0, 66.155),
            (2014.999, 67.619707),
            (2015.0, 67.62),
            (2017.0, 68.41134188381358),
            (2020.0, 69.37697312914538),
            (2023.0, 69.29761103397021),
            (2026.0, 69.0354672334697),
            (2026.5, 69.14),
            (2027.0, 69.28546094982505),
            (2030.0, 70.21847606910971),
            (2045.0, 76.43282247413141),
            (2100.0, 121.31021341515171),
            (3000.0, 5787.51292709445),
        ];

        for (year, expected) in cases {
            let actual = estimate(year).unwrap();
            assert!(
                (actual - expected).abs() < EPSILON,
                "year {year}: {actual} vs {expected}"
            );
        }
    }

    #[test]
    fn test_delta_t_recent_iers_observations() {
        // IERS 20u24 C04, downloaded 14 September 2026; ΔT = 32.184 + TAI-UTC - UT1-UTC.
        let cases = [
            (2015.0, 67.6439282),
            (2017.0, 68.5927130),
            (2017.4246575342465, 68.8085579),
            (2020.0, 69.3611665),
            (2023.0, 69.2038475),
            (2026.0, 69.1099131),
            (2026.4986301369863, 69.1691721),
        ];

        for (year, observed) in cases {
            let estimated = estimate(year).unwrap();
            assert!(
                (estimated - observed).abs() < 0.22,
                "year {year}: {estimated} vs {observed}"
            );
        }
    }

    #[test]
    fn test_delta_t_continuous_at_updated_boundaries() {
        for (year, expected) in [(2015.0_f64, 67.62), (2026.5, 69.14)] {
            for adjacent_year in [
                f64::from_bits(year.to_bits() - 1),
                year,
                f64::from_bits(year.to_bits() + 1),
            ] {
                let actual = estimate(adjacent_year).unwrap();
                assert!(
                    (actual - expected).abs() < 1e-12,
                    "year {adjacent_year}: {actual}"
                );
            }
        }
    }

    #[test]
    fn test_delta_t_smooth_joins_at_updated_boundaries() {
        for year in [2015.0, 2026.5] {
            let step = 1e-4;
            let value = estimate(year).unwrap();
            let left_slope = (value - estimate(year - step).unwrap()) / step;
            let right_slope = (estimate(year + step).unwrap() - value) / step;
            assert!(
                (left_slope - right_slope).abs() < 2e-5,
                "year {year}: left slope {left_slope} vs right slope {right_slope}"
            );
        }
    }

    #[test]
    fn test_delta_t_date_helpers_at_updated_branches() {
        for (year, month) in [(2014, 12), (2015, 1), (2026, 6), (2026, 7)] {
            let decimal_year = f64::from(year) + (f64::from(month) - 0.5) / 12.0;
            let expected = estimate(decimal_year).unwrap();
            assert_eq!(estimate_from_date(year, month).unwrap(), expected);

            #[cfg(feature = "chrono")]
            for day in [1, 28] {
                let date = chrono::NaiveDate::from_ymd_opt(year, month, day).unwrap();
                assert_eq!(estimate_from_date_like(date).unwrap(), expected);
            }
        }
    }

    #[test]
    fn test_delta_t_historical_estimates() {
        let delta_t_1900 = estimate(1900.0).unwrap();
        let delta_t_1950 = estimate(1950.0).unwrap();

        assert!(delta_t_1900 < 0.0); // Negative in early 20th century
        assert!(delta_t_1950 > 25.0 && delta_t_1950 < 35.0);
    }

    #[test]
    fn test_delta_t_boundary_conditions() {
        // Test edge cases
        assert!(estimate(-500.0).is_ok());
        assert!(estimate(3000.0).is_ok());
        assert!(estimate(-501.0).is_err());
        assert!(estimate(3001.0).is_err()); // Should fail beyond 3000
        assert!(estimate(f64::from_bits(3000.0_f64.to_bits() + 1)).is_err());
    }

    #[test]
    fn test_delta_t_from_date() {
        let delta_t = estimate_from_date(2024, 6).unwrap();
        let delta_t_decimal = estimate(2024.5 - 1.0 / 24.0).unwrap(); // June = month 6, so (6-0.5)/12 ≈ 0.458

        // Should be very close
        assert!((delta_t - delta_t_decimal).abs() < 0.01);

        // Test invalid month
        assert!(estimate_from_date(2024, 13).is_err());
        assert!(estimate_from_date(2024, 0).is_err());
    }

    #[test]
    #[cfg(feature = "chrono")]
    fn test_delta_t_from_date_like() {
        use chrono::{DateTime, FixedOffset, NaiveDate, Utc};

        // Test with DateTime<FixedOffset>
        let datetime_fixed = "2024-06-15T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let delta_t_fixed = estimate_from_date_like(datetime_fixed).unwrap();

        // Test with DateTime<Utc>
        let datetime_utc = "2024-06-15T19:00:00Z".parse::<DateTime<Utc>>().unwrap();
        let delta_t_utc = estimate_from_date_like(datetime_utc).unwrap();

        // Test with NaiveDate
        let naive_date = NaiveDate::from_ymd_opt(2024, 6, 15).unwrap();
        let delta_t_naive_date = estimate_from_date_like(naive_date).unwrap();

        // Test with NaiveDateTime
        let naive_datetime = naive_date.and_hms_opt(12, 0, 0).unwrap();
        let delta_t_naive_datetime = estimate_from_date_like(naive_datetime).unwrap();

        // Should all be identical since we only use year/month
        assert_eq!(delta_t_fixed, delta_t_utc);
        assert_eq!(delta_t_fixed, delta_t_naive_date);
        assert_eq!(delta_t_fixed, delta_t_naive_datetime);

        // Should match estimate_from_date
        let delta_t_date = estimate_from_date(2024, 6).unwrap();
        assert_eq!(delta_t_fixed, delta_t_date);

        // Verify reasonable range for 2024
        assert!(delta_t_fixed > 60.0 && delta_t_fixed < 80.0);
    }
}
