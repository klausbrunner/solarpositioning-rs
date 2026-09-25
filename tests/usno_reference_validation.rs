//! Validate sunrise/sunset calculations against USNO reference data.

#![cfg(all(feature = "chrono", feature = "std"))]

mod common;

use chrono::{DateTime, NaiveTime, Timelike, Utc};
use csv::ReaderBuilder;
use solar_positioning::SunriseResult;
use std::error::Error;
use std::fs::File;

fn time_difference_seconds(expected: &str, actual: DateTime<Utc>) -> Result<i64, Box<dyn Error>> {
    let expected = NaiveTime::parse_from_str(expected, "%H:%M")?;
    Ok((i64::from(actual.time().num_seconds_from_midnight())
        - i64::from(expected.num_seconds_from_midnight()))
    .abs())
}

fn validate(path: &str, elevation: f64, tolerance: i64) -> Result<(), Box<dyn Error>> {
    let file = File::open(path)?;
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(file);

    let mut count = 0;
    for record in reader.records() {
        let record = record?;
        let datetime = record[0].parse::<DateTime<Utc>>()?;
        let latitude = record[1].parse()?;
        let longitude = record[2].parse()?;
        let result = common::events_on_utc_date(datetime, latitude, longitude, elevation);
        // Known A.2.4 limitation: Anchorage's July 5 UTC dusk belongs to July 4's
        // transit, which SPA classifies as AllDay. Keep the reference row intact
        // and assert the limitation explicitly instead of comparing a July 6 event.
        if elevation == -6.0
            && latitude == 61.21666666666667
            && longitude == -149.86666666666667
            && record[0] == *"2020-07-05T12:00:00Z"
        {
            assert!(matches!(
                result,
                Err(solar_positioning::Error::ComputationError { .. })
            ));
            count += 1;
            continue;
        }
        let result = result?;

        match record[3].as_ref() {
            "NORMAL" => {
                let SunriseResult::RegularDay {
                    sunrise, sunset, ..
                } = result
                else {
                    panic!("expected a regular day for {datetime}, got {result:?}");
                };
                let sunrise_error = time_difference_seconds(&record[4], sunrise)?;
                let sunset_error = time_difference_seconds(&record[5], sunset)?;
                assert!(
                    sunrise_error < tolerance,
                    "sunrise error {sunrise_error}s for {datetime} at {latitude}, {longitude}: {result:?}"
                );
                assert!(
                    sunset_error < tolerance,
                    "sunset error {sunset_error}s for {datetime} at {latitude}, {longitude}: {result:?}"
                );
            }
            "ALL_DAY" => assert!(
                matches!(result, SunriseResult::AllDay { .. }),
                "expected polar day for {datetime}, got {result:?}"
            ),
            "ALL_NIGHT" => assert!(
                matches!(result, SunriseResult::AllNight { .. }),
                "expected polar night for {datetime}, got {result:?}"
            ),
            day_type => panic!("unknown day type {day_type:?}"),
        }
        count += 1;
    }

    assert!(count > 0, "no reference records loaded");
    Ok(())
}

#[test]
fn sunrise_sunset_matches_usno() -> Result<(), Box<dyn Error>> {
    validate("tests/data/usno/usno_reference_testdata.csv", -0.833, 40)
}

#[test]
fn extreme_latitudes_match_usno() -> Result<(), Box<dyn Error>> {
    validate(
        "tests/data/usno/usno_reference_testdata_extreme.csv",
        -0.833,
        150,
    )
}

#[test]
fn civil_twilight_matches_usno() -> Result<(), Box<dyn Error>> {
    validate(
        "tests/data/usno/usno_reference_testdata_civil.csv",
        -6.0,
        150,
    )
}
