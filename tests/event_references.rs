//! Independent reference events, shared with the Java implementation.
#![cfg(all(feature = "chrono", feature = "alloc"))]

use chrono::{DateTime, Duration, NaiveDate, NaiveTime, Utc};
use solar_positioning::SolarPositions;
use solar_positioning::{Horizon, Location, SolarEvents, time::JulianDate};

const DELTA_T: f64 = 69.184;

fn compare(expected: &[DateTime<Utc>], actual: &[DateTime<Utc>], tolerance: f64) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "expected {expected:?}, got {actual:?}"
    );
    for (reference, time) in expected.iter().zip(actual) {
        let difference = (*time - *reference).num_microseconds().unwrap().abs() as f64 / 1e6;
        assert!(
            difference <= tolerance,
            "{reference}: got {time}, error {difference}s > {tolerance}s"
        );
    }
}

fn parse_events(text: &str) -> Vec<DateTime<Utc>> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(';').map(|s| s.parse().unwrap()).collect()
    }
}

fn altitude(time: DateTime<Utc>, latitude: f64, longitude: f64) -> f64 {
    let jd = 2_440_587.5
        + time.timestamp() as f64 / 86400.0
        + f64::from(time.timestamp_subsec_nanos()) / 86400e9;
    SolarPositions::new()
        .at_from_julian(
            JulianDate::new(jd, DELTA_T).unwrap(),
            Location {
                latitude,
                longitude,
            },
            0.0,
            None,
        )
        .unwrap()
        .elevation_angle()
}

#[test]
fn matches_independent_jpl_events() {
    let calculator = SolarEvents::new();
    let mut reader = csv::Reader::from_path("tests/data/jpl_events.csv").unwrap();
    let mut count = 0;
    for row in reader.records() {
        let row = row.unwrap();
        let start: DateTime<Utc> = row[0].parse().unwrap();
        let end: DateTime<Utc> = row[1].parse().unwrap();
        let latitude = row[2].parse().unwrap();
        let longitude = row[3].parse().unwrap();
        let horizon = Horizon::Custom(row[4].parse().unwrap());
        let day = calculator
            .for_date(
                start.date_naive(),
                &Utc,
                Location {
                    latitude,
                    longitude,
                },
                DELTA_T,
                horizon,
            )
            .unwrap();
        assert_eq!(day.start, start);
        assert_eq!(day.end, end);
        for (field, actual) in [(&row[5], &day.rises), (&row[6], &day.sets)] {
            let expected = parse_events(field);
            // SPA's uncertainty is angular. At shallow crossings, convert it to
            // a time tolerance using local elevation speed, without relaxing counts.
            let tolerance = expected.iter().fold(0.002_f64, |limit, &time| {
                let speed = (altitude(time + Duration::seconds(1), latitude, longitude)
                    - altitude(time - Duration::seconds(1), latitude, longitude))
                .abs()
                    / 2.0;
                limit.max(0.0003 / speed + 0.002)
            });
            compare(&expected, actual, tolerance);
        }
        compare(&parse_events(&row[7]), &day.transits, 1.0);
        count += 1;
    }
    assert_eq!(count, 56);
}

fn compare_minute(date: NaiveDate, expected: &str, actual: &[DateTime<Utc>]) {
    let expected: Vec<_> = if expected.is_empty() {
        Vec::new()
    } else {
        vec![
            date.and_time(NaiveTime::parse_from_str(expected, "%H:%M").unwrap())
                .and_utc(),
        ]
    };
    compare(&expected, actual, 90.0);
}

#[test]
fn matches_usno_sunrise_and_twilight_tables() {
    let calculator = SolarEvents::new();
    for (file, horizon) in [
        ("usno_reference_testdata.csv", Horizon::SunriseSunset),
        (
            "usno_reference_testdata_extreme.csv",
            Horizon::SunriseSunset,
        ),
        ("usno_reference_testdata_civil.csv", Horizon::CivilTwilight),
    ] {
        let mut reader = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .has_headers(false)
            .from_path(format!("tests/data/usno/{file}"))
            .unwrap();
        let mut count = 0;
        for row in reader.records() {
            let row = row.unwrap();
            let date = row[0].parse::<DateTime<Utc>>().unwrap().date_naive();
            let day = calculator
                .for_date(
                    date,
                    &Utc,
                    Location {
                        latitude: row[1].parse().unwrap(),
                        longitude: row[2].parse().unwrap(),
                    },
                    DELTA_T,
                    horizon,
                )
                .unwrap();
            compare_minute(date, &row[4], &day.rises);
            compare_minute(date, &row[5], &day.sets);
            assert_eq!(day.always_above(), &row[3] == "ALL_DAY", "{row:?}");
            assert_eq!(day.always_below(), &row[3] == "ALL_NIGHT", "{row:?}");
            count += 1;
        }
        assert!(count > 0, "empty fixture {file}");
    }
}
