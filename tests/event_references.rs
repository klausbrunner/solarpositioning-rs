//! Independent reference events, shared with the Java implementation.
#![cfg(all(feature = "jiff", feature = "alloc"))]

use jiff::{
    SignedDuration, Timestamp, ZonedDateTime,
    civil::{Date, Time},
    tz::TimeZone,
};
use solar_positioning::{Horizon, Location, SolarEvents, SolarPositions};

const DELTA_T: f64 = 69.184;

fn compare(expected: &[Timestamp], actual: &[ZonedDateTime], tolerance: f64) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "expected {expected:?}, got {actual:?}"
    );
    for (reference, time) in expected.iter().zip(actual) {
        let difference = (time.timestamp() - *reference).as_secs_f64().abs();
        assert!(
            difference <= tolerance,
            "{reference}: got {time}, error {difference}s > {tolerance}s"
        );
    }
}

fn parse_events(text: &str) -> Vec<Timestamp> {
    if text.is_empty() {
        Vec::new()
    } else {
        text.split(';').map(|s| s.parse().unwrap()).collect()
    }
}

fn altitude(time: Timestamp, location: Location) -> f64 {
    SolarPositions::new()
        .at(time, location, 0.0, DELTA_T, None)
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
        let start: Timestamp = row[0].parse().unwrap();
        let end: Timestamp = row[1].parse().unwrap();
        let location = Location {
            latitude: row[2].parse().unwrap(),
            longitude: row[3].parse().unwrap(),
        };
        let horizon = Horizon::Custom(row[4].parse().unwrap());
        let day = calculator
            .for_date(
                TimeZone::UTC.to_datetime(start).date(),
                &TimeZone::UTC,
                location,
                DELTA_T,
                horizon,
            )
            .unwrap();
        assert_eq!(day.start.timestamp(), start);
        assert_eq!(day.end.timestamp(), end);
        for (field, actual) in [(&row[5], &day.rises), (&row[6], &day.sets)] {
            let expected = parse_events(field);
            // SPA's uncertainty is angular. At shallow crossings, convert it to
            // a time tolerance using local elevation speed, without relaxing counts.
            let tolerance = expected.iter().fold(0.101_f64, |limit, &time| {
                let speed = (altitude(time + SignedDuration::from_secs(1), location)
                    - altitude(time - SignedDuration::from_secs(1), location))
                .abs()
                    / 2.0;
                limit.max(0.0003 / speed + 0.101)
            });
            compare(&expected, actual, tolerance);
        }
        compare(&parse_events(&row[7]), &day.transits, 1.0);
        count += 1;
    }
    assert_eq!(count, 56);
}

fn compare_minute(date: Date, expected: &str, actual: &[ZonedDateTime]) {
    let expected = (!expected.is_empty()).then(|| {
        TimeZone::UTC
            .to_timestamp(date.to_datetime(expected.parse::<Time>().unwrap()))
            .unwrap()
    });
    compare(expected.as_slice(), actual, 90.0);
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
            let date = TimeZone::UTC.to_datetime(row[0].parse().unwrap()).date();
            let day = calculator
                .for_date(
                    date,
                    &TimeZone::UTC,
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
