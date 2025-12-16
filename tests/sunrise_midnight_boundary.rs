use chrono::{DateTime, FixedOffset, NaiveDate};
use solar_positioning::{spa, time::DeltaT, Horizon, SunriseResult};

// When sunrise/sunset events occur near local midnight, a naive "same-local-date" mapping can
// skip/duplicate events. The chrono API uses the local calendar date as the SPA UTC day and
// corrects sunrise when it ends up after transit.
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
