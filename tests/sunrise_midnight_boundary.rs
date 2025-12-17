use chrono::{DateTime, FixedOffset, NaiveDate};
use solar_positioning::{spa, time::DeltaT, Horizon, SunriseResult};

// When sunrise/sunset events occur near local midnight, a naive UTC-day mapping can
// skip/duplicate events. The chrono API picks the internal UTC day so that transit lands on the
// requested local calendar date, and corrects sunrise when it ends up after transit.
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
