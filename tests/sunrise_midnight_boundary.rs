#![cfg(feature = "chrono")]

use chrono::{DateTime, Duration, FixedOffset, NaiveDate, TimeZone};
use solar_positioning::{spa, time::DeltaT, Horizon, SunriseResult};

// When sunrise/sunset events occur near local midnight, a naive UTC-day mapping can
// skip/duplicate events. The chrono API selects a transit near local noon, retaining
// event day offsets throughout the calculation.
#[test]
fn sunrise_near_local_midnight_is_not_skipped() {
    let latitude = 49.60139790853522;
    let longitude = 171.01752655220554;
    let horizon = Horizon::AstronomicalTwilight;

    for day in [1_u32, 2, 3] {
        let date = format!("1986-06-0{day}T00:00:00+11:00")
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let delta_t = DeltaT::estimate_from_date_like(date).unwrap();

        let result =
            spa::sunrise_sunset_for_horizon(date, latitude, longitude, delta_t, horizon).unwrap();

        let sunrise = match result {
            SunriseResult::RegularDay { sunrise, .. } => sunrise,
            _ => panic!("expected regular day"),
        };

        let minutes_from_local_midnight = sunrise.signed_duration_since(date).num_minutes();
        assert!(
            (-60..=60).contains(&minutes_from_local_midnight),
            "day={day} sunrise={sunrise} minutes_from_midnight={minutes_from_local_midnight}"
        );

        if day == 3 {
            assert!(
                minutes_from_local_midnight < 0,
                "expected day 3 sunrise to occur before local midnight, got {sunrise}"
            );
            assert_eq!(
                sunrise.date_naive(),
                NaiveDate::from_ymd_opt(1986, 6, 2).unwrap()
            );
        }
    }
}

#[test]
fn sunrise_sunset_multiple_matches_midnight_handling() {
    let latitude = 49.60139790853522;
    let longitude = 171.01752655220554;
    let horizon = Horizon::AstronomicalTwilight;

    let date = "1986-06-03T00:00:00+11:00"
        .parse::<DateTime<FixedOffset>>()
        .unwrap();
    let delta_t = DeltaT::estimate_from_date_like(date).unwrap();

    let results: Vec<_> =
        spa::sunrise_sunset_multiple(date, latitude, longitude, delta_t, [horizon])
            .collect::<solar_positioning::Result<Vec<_>>>()
            .unwrap();

    assert_eq!(results.len(), 1);
    let (_h, result) = &results[0];

    let sunrise = match result {
        SunriseResult::RegularDay { sunrise, .. } => *sunrise,
        _ => panic!("expected regular day"),
    };

    let minutes_from_local_midnight = sunrise.signed_duration_since(date).num_minutes();
    assert!(
        (-60..=60).contains(&minutes_from_local_midnight),
        "sunrise={sunrise} minutes_from_midnight={minutes_from_local_midnight}"
    );
    assert!(minutes_from_local_midnight < 0);
    assert_eq!(
        sunrise.date_naive(),
        NaiveDate::from_ymd_opt(1986, 6, 2).unwrap()
    );
}

#[test]
fn sunrise_near_antimeridian_is_not_shifted_to_next_day() {
    // Regression for timezone offsets near the antimeridian where sunrise can fall late in UTC.
    let latitude = -9.459488331200241;
    let longitude = 177.60664224032377;

    let date = "1997-11-01T00:00:00+12:00"
        .parse::<DateTime<FixedOffset>>()
        .unwrap();
    let delta_t = DeltaT::estimate_from_date_like(date).unwrap();

    let result =
        spa::sunrise_sunset_for_horizon(date, latitude, longitude, delta_t, Horizon::SunriseSunset)
            .unwrap();

    let (sunrise, transit) = match result {
        SunriseResult::RegularDay {
            sunrise, transit, ..
        } => (sunrise, transit),
        _ => panic!("expected regular day"),
    };

    assert!(sunrise < transit, "sunrise={sunrise} transit={transit}");
    assert_eq!(
        sunrise.date_naive(),
        NaiveDate::from_ymd_opt(1997, 11, 1).unwrap()
    );
}

#[test]
fn sunrise_stays_on_local_date_for_plus_five_offset() {
    let latitude = -29.807961253888443;
    let longitude = 80.510704531778;

    let date = "2008-10-16T00:00:00+05:00"
        .parse::<DateTime<FixedOffset>>()
        .unwrap();
    let delta_t = DeltaT::estimate_from_date_like(date).unwrap();

    let result =
        spa::sunrise_sunset_for_horizon(date, latitude, longitude, delta_t, Horizon::SunriseSunset)
            .unwrap();

    let (sunrise, transit, sunset) = match result {
        SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        } => (sunrise, transit, sunset),
        _ => panic!("expected regular day"),
    };

    let local_date = date.date_naive();
    assert_eq!(
        transit.date_naive(),
        local_date,
        "transit={transit} sunrise={sunrise} sunset={sunset}"
    );
    assert_eq!(
        sunrise.date_naive(),
        local_date,
        "transit={transit} sunrise={sunrise} sunset={sunset}"
    );
    assert!(sunrise < transit, "sunrise={sunrise} transit={transit}");
    assert!(transit < sunset, "transit={transit} sunset={sunset}");
}

#[test]
fn transit_stays_on_requested_day_across_antimeridian() {
    // Around the equation-of-time zero, wrapped transit jumps between UTC days.
    for (longitude, offset) in [(-180.0, "-12:00"), (180.0, "+12:00")] {
        for day in 12..=18 {
            let date = format!("2024-04-{day:02}T00:00:00{offset}")
                .parse::<DateTime<FixedOffset>>()
                .unwrap();
            for latitude in [-60.0, 0.0, 60.0] {
                for (horizon, result) in spa::sunrise_sunset_multiple(
                    date,
                    latitude,
                    longitude,
                    69.184,
                    [Horizon::SunriseSunset, Horizon::CivilTwilight],
                )
                .map(Result::unwrap)
                {
                    assert_eq!(result.transit().date_naive(), date.date_naive());
                    assert_eq!(result, spa::sunrise_sunset_for_horizon(
                        date, latitude, longitude, 69.184, horizon,
                    ).unwrap());
                    let SunriseResult::RegularDay {
                        sunrise,
                        transit,
                        sunset,
                    } = result
                    else {
                        panic!("expected regular day");
                    };
                    assert!(sunrise < transit && transit < sunset);
                    assert!((transit - sunrise).num_hours() < 24);
                    assert!((sunset - transit).num_hours() < 24);
                }
            }
        }
    }
}

#[test]
fn sunset_uses_its_actual_day_across_utc_midnight() {
    // Skyfield 1.55 / JPL DE440s: sea level, delta_t=69.184s,
    // geometric solar-centre horizon -0.83337 degrees. Previously 9m16s late.
    let date = "2024-03-20T00:00:00-08:00"
        .parse::<DateTime<FixedOffset>>()
        .unwrap();
    let expected = "2024-03-21T02:17:36Z"
        .parse::<DateTime<FixedOffset>>()
        .unwrap();
    let result =
        spa::sunrise_sunset_for_horizon(date, -80.0, -120.0, 69.184, Horizon::SunriseSunset)
            .unwrap();
    assert!((*result.sunset().unwrap() - expected).num_seconds().abs() < 5);
}

// Eight civil dates without a transit, then two with two transits. The expected
// SPA candidate nearest local noon defines the cycle, not an event shifted by 24h.
#[test]
fn transit_selection_handles_missing_and_double_transits() {
    for (date, longitude, expected) in [
        ("2020-06-10", 179.9, "2020-06-11T00:00:03.291Z"),
        ("2020-12-24", 179.9, "2020-12-23T23:59:56.717Z"),
        ("2020-06-14", -179.9, "2020-06-15T00:00:05.462Z"),
        ("2020-12-25", -179.9, "2020-12-26T00:00:08.100Z"),
        ("2025-06-10", 179.9, "2025-06-11T00:00:00.776Z"),
        ("2025-12-24", 179.9, "2025-12-23T23:59:49.820Z"),
        ("2025-06-14", -179.9, "2025-06-15T00:00:02.694Z"),
        ("2025-12-25", -179.9, "2025-12-26T00:00:01.383Z"),
        ("2020-04-16", 179.9, "2020-04-16T00:00:12.776Z"),
        ("2020-09-02", 179.9, "2020-09-02T23:59:45.609Z"),
    ] {
        let day = format!("{date}T00:00:00Z")
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let expected = expected.parse::<DateTime<FixedOffset>>().unwrap();
        for offset in -1..=1 {
            let query = day + Duration::days(offset);
            let delta_t = DeltaT::estimate_from_date_like(query).unwrap();
            for (horizon, result) in spa::sunrise_sunset_multiple(
                query,
                0.0,
                longitude,
                delta_t,
                [Horizon::SunriseSunset, Horizon::CivilTwilight],
            )
            .map(Result::unwrap)
            {
                if offset == 0 {
                    assert!(
                        (*result.transit() - expected).num_milliseconds().abs() <= 1,
                        "{query}, {longitude}: {:?}",
                        result.transit()
                    );
                } else {
                    assert_eq!(result.transit().date_naive(), query.date_naive());
                }
                assert_eq!(
                    result,
                    spa::sunrise_sunset_for_horizon(query, 0.0, longitude, delta_t, horizon)
                        .unwrap()
                );
                let SunriseResult::RegularDay {
                    sunrise,
                    transit,
                    sunset,
                } = result
                else {
                    panic!("expected regular solar cycle");
                };
                assert!(sunrise < transit && transit < sunset);
            }
        }
    }
}

#[test]
fn transit_selection_uses_local_clock_across_dst() {
    for (month, day, offset) in [(3, 31, 7200), (10, 27, 3600)] {
        let query = chrono_tz::Europe::Berlin
            .with_ymd_and_hms(2024, month, day, 0, 0, 0)
            .unwrap();
        let result =
            spa::sunrise_sunset_for_horizon(query, 52.0, 13.4, 69.184, Horizon::SunriseSunset)
                .unwrap();
        assert_eq!(result.transit().date_naive(), query.date_naive());
        assert_eq!(
            result.transit().fixed_offset().offset().local_minus_utc(),
            offset
        );
        assert_eq!(
            result,
            spa::sunrise_sunset_for_horizon(
                query + Duration::hours(12),
                52.0,
                13.4,
                69.184,
                Horizon::SunriseSunset
            )
            .unwrap()
        );
    }
}

#[test]
fn transit_selection_handles_clock_rollbacks_and_calendar_gaps() {
    use chrono_tz::{America::Adak, UTC};
    for (date, zone, longitude, expected) in [
        ("1867-10-18", Adak, 0.0, "1867-10-19T11:45:05.512Z"),
        ("1867-10-19", Adak, -176.64, "1867-10-18T23:31:44.755Z"),
        ("1582-10-04", UTC, 0.0, "1582-10-04T11:46:11.877Z"),
        ("1582-10-15", UTC, -180.0, "1582-10-15T23:45:52.632Z"),
    ] {
        let local = date
            .parse::<NaiveDate>()
            .unwrap()
            .and_hms_opt(12, 0, 0)
            .unwrap();
        let query = zone.from_local_datetime(&local).earliest().unwrap();
        let expected = expected.parse::<DateTime<FixedOffset>>().unwrap();
        for (horizon, result) in spa::sunrise_sunset_multiple(
            query,
            0.0,
            longitude,
            69.184,
            [Horizon::SunriseSunset, Horizon::CivilTwilight],
        )
        .map(Result::unwrap)
        {
            assert!(
                (result.transit().timestamp_millis() - expected.timestamp_millis()).abs() <= 1,
                "{date}, {zone}, {longitude}: {:?}",
                result.transit()
            );
            assert_eq!(
                result,
                spa::sunrise_sunset_for_horizon(query, 0.0, longitude, 69.184, horizon).unwrap()
            );
        }
    }
    // Skipping invalid neighbours must not accept an invalid requested date.
    let invalid = "1582-10-10T12:00:00Z"
        .parse::<DateTime<FixedOffset>>()
        .unwrap();
    assert!(matches!(
        spa::sunrise_sunset_for_horizon(invalid, 0.0, 0.0, 69.184, Horizon::SunriseSunset),
        Err(solar_positioning::Error::InvalidDateTime { .. })
    ));
}
