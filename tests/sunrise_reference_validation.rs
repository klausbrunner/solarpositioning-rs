//! Validate sunrise/sunset calculations against SPA reference data.

#![cfg(all(feature = "chrono", feature = "std"))]

mod common;

use chrono::{DateTime, NaiveTime, Timelike, Utc};
use csv::ReaderBuilder;
use solar_positioning::{spa, SunriseResult};
use std::error::Error;
use std::fs::File;

const TRANSIT_TOLERANCE_SECONDS: i64 = 1;
const HORIZON_TOLERANCE_SECONDS: i64 = 120;

fn time_difference_seconds(expected: &str, actual: DateTime<Utc>) -> Result<i64, Box<dyn Error>> {
    let expected = NaiveTime::parse_from_str(expected, "%H:%M:%S")?;
    Ok((i64::from(actual.time().num_seconds_from_midnight())
        - i64::from(expected.num_seconds_from_midnight()))
    .abs())
}

#[test]
fn sunrise_sunset_matches_spa_reference_data() -> Result<(), Box<dyn Error>> {
    let file = File::open("tests/data/test/sunrise/spa_reference_testdata.csv")?;
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

        match common::events_on_utc_date(datetime, latitude, longitude, -0.833)? {
            SunriseResult::RegularDay {
                sunrise,
                transit,
                sunset,
            } => {
                let sunrise_error = time_difference_seconds(&record[3], sunrise)?;
                let transit_error = time_difference_seconds(&record[4], transit)?;
                let sunset_error = time_difference_seconds(&record[5], sunset)?;

                assert!(
                    sunrise_error < HORIZON_TOLERANCE_SECONDS,
                    "sunrise error {sunrise_error}s for {datetime}"
                );
                assert!(
                    transit_error < TRANSIT_TOLERANCE_SECONDS,
                    "transit error {transit_error}s for {datetime}"
                );
                assert!(
                    sunset_error < HORIZON_TOLERANCE_SECONDS,
                    "sunset error {sunset_error}s for {datetime}"
                );
            }
            SunriseResult::AllDay { .. } | SunriseResult::AllNight { .. } => {
                assert!(
                    record.iter().skip(3).all(str::is_empty),
                    "unexpected polar day or night for {datetime}"
                );
            }
        }
        count += 1;
    }

    assert!(count > 0, "no reference records loaded");
    Ok(())
}

#[test]
fn polar_transit_matches_utc_api() -> Result<(), Box<dyn Error>> {
    let latitude = 78.0;
    let longitude = 15.0;
    let date = "2024-01-01T12:00:00+01:00".parse::<DateTime<chrono::FixedOffset>>()?;

    let SunriseResult::AllNight { transit } =
        spa::sunrise_sunset(date, latitude, longitude, 0.0, -0.833)?
    else {
        panic!("expected polar night");
    };
    let SunriseResult::AllNight { transit: expected } =
        spa::sunrise_sunset_utc(2024, 1, 1, latitude, longitude, 0.0, -0.833)?
    else {
        panic!("expected polar night from UTC API");
    };

    let expected = "2024-01-01T00:00:00Z".parse::<DateTime<Utc>>()?
        + chrono::Duration::milliseconds((expected.hours() * 3_600_000.0) as i64);
    assert_eq!(transit.with_timezone(&Utc), expected);

    let position = spa::solar_position(transit, latitude, longitude, 0.0, 0.0, None)?;
    assert!((position.azimuth() - 180.0).abs() < 0.001);
    Ok(())
}
