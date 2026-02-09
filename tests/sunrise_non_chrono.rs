//! Tests for non-chrono sunrise/sunset API

use solar_positioning::{spa, Horizon, HoursUtc, SunriseResult};

#[test]
fn test_sunrise_sunset_utc_basic() {
    // San Francisco, June 21, 2023
    let result = spa::sunrise_sunset_utc(2023, 6, 21, 37.7749, -122.4194, 69.0, -0.833).unwrap();

    if let SunriseResult::RegularDay {
        sunrise,
        transit,
        sunset,
    } = result
    {
        // Sunrise should be around 12:48 UTC (05:48 PDT)
        assert!((sunrise.hours() - 12.8).abs() < 0.5);
        // Transit should be around 20:11 UTC (13:11 PDT)
        assert!((transit.hours() - 20.2).abs() < 0.5);
        // Sunset should be around 03:35 UTC next day (20:35 PDT same day)
        assert!((sunset.hours() - 3.58).abs() < 0.5);
    } else {
        panic!("Expected RegularDay result");
    }
}

#[test]
fn test_sunrise_sunset_utc_for_horizon() {
    // Test with different horizon types
    let horizons = [
        Horizon::SunriseSunset,
        Horizon::CivilTwilight,
        Horizon::NauticalTwilight,
        Horizon::AstronomicalTwilight,
    ];

    for horizon in horizons {
        let result =
            spa::sunrise_sunset_utc_for_horizon(2023, 6, 21, 37.7749, -122.4194, 69.0, horizon)
                .unwrap();

        // All should return RegularDay for San Francisco in June
        assert!(matches!(result, SunriseResult::RegularDay { .. }));
    }
}

#[test]
fn test_hours_utc_day_and_hours() {
    // Test current day
    let h1 = HoursUtc::from_hours(12.5);
    let (day, hours) = h1.day_and_hours();
    assert_eq!(day, 0);
    assert!((hours - 12.5).abs() < 1e-10);

    // Test next day
    let h2 = HoursUtc::from_hours(25.5);
    let (day, hours) = h2.day_and_hours();
    assert_eq!(day, 1);
    assert!((hours - 1.5).abs() < 1e-10);

    // Test previous day
    let h3 = HoursUtc::from_hours(-0.5);
    let (day, hours) = h3.day_and_hours();
    assert_eq!(day, -1);
    assert!((hours - 23.5).abs() < 1e-10);

    // Multi-day positive offset
    let h4 = HoursUtc::from_hours(49.25);
    let (day, hours) = h4.day_and_hours();
    assert_eq!(day, 2);
    assert!((hours - 1.25).abs() < 1e-10);

    // Multi-day negative offset
    let h5 = HoursUtc::from_hours(-50.5);
    let (day, hours) = h5.day_and_hours();
    assert_eq!(day, -3);
    assert!((hours - 21.5).abs() < 1e-10);
}

#[test]
fn test_polar_regions() {
    // Test polar day (Svalbard in summer)
    let result = spa::sunrise_sunset_utc_for_horizon(
        2023,
        6,
        21,
        78.0, // Svalbard latitude
        15.0,
        69.0,
        Horizon::SunriseSunset,
    )
    .unwrap();

    assert!(matches!(result, SunriseResult::AllDay { .. }));

    // Test polar night (Svalbard in winter)
    let result = spa::sunrise_sunset_utc_for_horizon(
        2023,
        12,
        21,
        78.0, // Svalbard latitude
        15.0,
        69.0,
        Horizon::SunriseSunset,
    )
    .unwrap();

    assert!(matches!(result, SunriseResult::AllNight { .. }));
}

#[test]
fn test_invalid_inputs() {
    // Invalid latitude
    let result = spa::sunrise_sunset_utc(2023, 6, 21, 91.0, 0.0, 69.0, -0.833);
    assert!(result.is_err());

    // Invalid longitude
    let result = spa::sunrise_sunset_utc(2023, 6, 21, 0.0, 181.0, 69.0, -0.833);
    assert!(result.is_err());

    // Invalid month
    let result = spa::sunrise_sunset_utc(2023, 13, 21, 0.0, 0.0, 69.0, -0.833);
    assert!(result.is_err());

    // Invalid day
    let result = spa::sunrise_sunset_utc(2023, 6, 32, 0.0, 0.0, 69.0, -0.833);
    assert!(result.is_err());

    // Invalid delta_t (non-finite)
    let result = spa::sunrise_sunset_utc(2023, 6, 21, 0.0, 0.0, f64::NAN, -0.833);
    assert!(result.is_err());
    let result = spa::sunrise_sunset_utc(2023, 6, 21, 0.0, 0.0, f64::INFINITY, -0.833);
    assert!(result.is_err());
}

#[test]
fn test_consistency_across_api() {
    // Non-chrono API
    let result_utc = spa::sunrise_sunset_utc(2023, 6, 21, 40.0, -75.0, 69.0, -0.833).unwrap();

    // Should return RegularDay
    if let SunriseResult::RegularDay {
        sunrise: sunrise_utc,
        transit: transit_utc,
        sunset: sunset_utc,
    } = result_utc
    {
        // Check that times are reasonable
        assert!(sunrise_utc.hours() > 5.0 && sunrise_utc.hours() < 15.0);
        assert!(transit_utc.hours() > 10.0 && transit_utc.hours() < 20.0);
        assert!(sunset_utc.hours() >= 0.0 && sunset_utc.hours() < 10.0);
    } else {
        panic!("Expected RegularDay result");
    }
}

#[cfg(feature = "chrono")]
#[test]
fn test_chrono_vs_non_chrono_consistency() {
    use chrono::{FixedOffset, NaiveDate, TimeZone, Timelike};

    // Non-chrono API
    let result_utc =
        spa::sunrise_sunset_utc(2023, 6, 21, 37.7749, -122.4194, 69.0, -0.833).unwrap();

    // Chrono API
    let date = FixedOffset::east_opt(0)
        .unwrap()
        .from_local_datetime(
            &NaiveDate::from_ymd_opt(2023, 6, 21)
                .unwrap()
                .and_hms_opt(0, 0, 0)
                .unwrap(),
        )
        .unwrap();
    let result_chrono = spa::sunrise_sunset(date, 37.7749, -122.4194, 69.0, -0.833).unwrap();

    // Both should return RegularDay
    match (result_utc, result_chrono) {
        (
            SunriseResult::RegularDay {
                sunrise: sunrise_utc,
                transit: transit_utc,
                sunset: sunset_utc,
            },
            SunriseResult::RegularDay {
                sunrise: sunrise_chrono,
                transit: transit_chrono,
                sunset: sunset_chrono,
            },
        ) => {
            // Convert chrono DateTime to hours
            let sunrise_hours = sunrise_chrono.hour() as f64
                + sunrise_chrono.minute() as f64 / 60.0
                + sunrise_chrono.second() as f64 / 3600.0;
            let transit_hours = transit_chrono.hour() as f64
                + transit_chrono.minute() as f64 / 60.0
                + transit_chrono.second() as f64 / 3600.0;
            let sunset_hours = sunset_chrono.hour() as f64
                + sunset_chrono.minute() as f64 / 60.0
                + sunset_chrono.second() as f64 / 3600.0;

            // Should match within 1 second (1/3600 hour)
            assert!((sunrise_utc.hours() - sunrise_hours).abs() < 1.0 / 3600.0);
            assert!((transit_utc.hours() - transit_hours).abs() < 1.0 / 3600.0);
            assert!((sunset_utc.hours() - sunset_hours).abs() < 1.0 / 3600.0);
        }
        _ => panic!("Expected RegularDay from both APIs"),
    }
}
