#![cfg(feature = "chrono")]

//! Regression tests covering tricky sunrise/sunset edge cases.

use chrono::{offset::LocalResult, NaiveDate, TimeZone};
use chrono_tz::TZ_VARIANTS;
use solar_positioning::{spa, Horizon, SunriseResult};

fn find_date_with_missing_midnight() -> Option<chrono::DateTime<chrono_tz::Tz>> {
    for tz in TZ_VARIANTS {
        for year in 2000..=2030 {
            for month in 1..=12 {
                for day in 1..=31 {
                    let Some(date) = NaiveDate::from_ymd_opt(year, month, day) else {
                        continue;
                    };

                    let midnight_local =
                        tz.from_local_datetime(&date.and_hms_opt(0, 0, 0).unwrap());
                    if !matches!(midnight_local, LocalResult::Single(_)) {
                        if let LocalResult::Single(midday) =
                            tz.from_local_datetime(&date.and_hms_opt(12, 0, 0).unwrap())
                        {
                            return Some(midday);
                        }
                    }
                }
            }
        }
    }
    None
}

#[test]
fn sunrise_handles_dates_without_local_midnight() {
    let datetime = find_date_with_missing_midnight()
        .expect("expected to find a timezone/day combination without midnight");

    // Equatorial coordinates keep behaviour simple and avoid polar day/night.
    let result = std::panic::catch_unwind(|| {
        spa::sunrise_sunset(
            datetime,
            0.0,
            0.0,
            0.0,
            Horizon::SunriseSunset.elevation_angle(),
        )
    });

    assert!(
        result.is_ok(),
        "sunrise calculation should not panic for {:?}",
        datetime
    );

    let sunrise_result = result.unwrap().expect("calculation should succeed");
    assert!(
        matches!(sunrise_result, SunriseResult::RegularDay { .. }),
        "expected regular sunrise/sunset result for equatorial coordinates"
    );
}

#[test]
fn sunrise_results_are_finite_near_polar_boundary() {
    let mut problematic_latitude = None;

    for latitude in (0..=480).map(|i| 65.0 + f64::from(i) * 0.05) {
        let outcome = spa::sunrise_sunset_utc(
            2023,
            6,
            21,
            latitude,
            0.0,
            69.0,
            Horizon::SunriseSunset.elevation_angle(),
        );

        if let Ok(SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        }) = outcome
        {
            if !sunrise.hours().is_finite()
                || !transit.hours().is_finite()
                || !sunset.hours().is_finite()
            {
                problematic_latitude = Some(latitude);
                break;
            }
        }
    }

    assert!(
        problematic_latitude.is_none(),
        "found non-finite sunrise result near latitude {:?}",
        problematic_latitude
    );
}
