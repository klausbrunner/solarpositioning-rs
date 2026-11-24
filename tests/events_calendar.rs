#![cfg(all(feature = "jiff", feature = "alloc"))]

use jiff::{
    SignedDuration, Timestamp,
    civil::{Date, date},
    tz::{TimeZone, offset},
};
use solar_positioning::SolarPositions;
use solar_positioning::{EventPosition, Horizon, Location, SolarEvents, time::JulianDate};

const DELTA_T: f64 = 69.184;
const HORIZONS: [Horizon; 4] = [
    Horizon::SunriseSunset,
    Horizon::CivilTwilight,
    Horizon::NauticalTwilight,
    Horizon::AstronomicalTwilight,
];

const ORIGIN: Location = Location {
    latitude: 0.0,
    longitude: 0.0,
};

fn instant(s: &str) -> Timestamp {
    s.parse().unwrap()
}
fn jd(time: Timestamp) -> f64 {
    2_440_587.5 + time.as_second() as f64 / 86400.0 + f64::from(time.subsec_nanosecond()) / 86400e9
}
fn compare(expected: &[Timestamp], actual: &[Timestamp]) {
    assert_eq!(
        expected.len(),
        actual.len(),
        "expected {expected:?}, got {actual:?}"
    );
    for (a, b) in expected.iter().zip(actual) {
        assert!(
            (*a - *b).abs() <= SignedDuration::from_millis(100),
            "{a} vs {b}"
        );
    }
}
fn find_all(
    mut start: Timestamp,
    end: Timestamp,
    next: impl Fn(Timestamp, Timestamp) -> solar_positioning::Result<Option<Timestamp>>,
) -> Vec<Timestamp> {
    let mut result = Vec::new();
    while let Some(time) = next(start, end).unwrap() {
        assert!(time > start && time <= end);
        result.push(time);
        assert!(result.len() < 10, "repeated event");
        start = time;
    }
    result
}

#[test]
fn handles_short_long_repeated_skipped_and_midnight_gap_dates() {
    for (date, zone, hours, count) in [
        ("2024-03-31", "Europe/Berlin", 23, 1),
        ("2024-10-27", "Europe/Berlin", 25, 1),
        ("1892-07-04", "Pacific/Apia", 48, 2),
        ("2011-12-30", "Pacific/Apia", 0, 0),
        ("2018-11-04", "America/Sao_Paulo", 23, 1),
    ] {
        let date = date.parse::<Date>().unwrap();
        let zone = TimeZone::get(zone).unwrap();
        let day = SolarEvents::new()
            .for_date(date, &zone, ORIGIN, DELTA_T, Horizon::SunriseSunset)
            .unwrap();
        assert_eq!(
            day.end.timestamp() - day.start.timestamp(),
            SignedDuration::from_hours(hours)
        );
        assert!(!day.always_above() && !day.always_below());
        for events in [&day.rises, &day.sets, &day.transits] {
            assert_eq!(events.len(), count, "{date} {zone:?}: {events:?}");
            assert!(events.windows(2).all(|pair| pair[0] < pair[1]));
            for time in events {
                assert_eq!(time.time_zone(), &zone);
                assert_eq!(time.date(), date);
                assert!(*time >= day.start && *time < day.end);
            }
        }
    }
}

#[test]
fn begins_at_the_first_valid_instant_when_a_gap_straddles_midnight() {
    // Local clocks jump from 23:30 to 00:30.
    let zone = TimeZone::posix("STD0DST,M3.2.0/-0:30,M11.1.0/2").unwrap();
    let date = date(2024, 3, 10);
    let day = SolarEvents::new()
        .for_date(date, &zone, ORIGIN, DELTA_T, Horizon::SunriseSunset)
        .unwrap();
    assert_eq!(day.start.datetime(), date.at(0, 30, 0, 0));
}

#[test]
fn horizons_share_transit_and_keep_independent_crossings() {
    let calculator = SolarEvents::new();
    let location = Location {
        latitude: 70.0,
        longitude: 0.0,
    };
    let date = date(2024, 12, 21);
    let horizons = [
        Horizon::SunriseSunset,
        Horizon::CivilTwilight,
        Horizon::Custom(-4.5),
        Horizon::CivilTwilight,
    ];
    let days = calculator
        .for_date_multiple(date, &TimeZone::UTC, location, DELTA_T, horizons)
        .unwrap();
    assert_eq!(days.len(), horizons.len());
    assert!(days[0].1.always_below());
    assert!(!days[1].1.always_below());
    assert_eq!(days[1].1.rises.len(), 1);
    for ((horizon, day), expected_horizon) in days.iter().zip(horizons) {
        assert_eq!(*horizon, expected_horizon);
        assert_eq!(
            *day,
            calculator
                .for_date(date, &TimeZone::UTC, location, DELTA_T, *horizon)
                .unwrap()
        );
        assert_eq!(day.transits, days[0].1.transits);
    }
}

#[test]
fn grena_events_follow_grena_positions_at_all_horizons() {
    let calculator = SolarEvents::grena3();
    for (date, lat, lon, rises, sets, transits) in [
        ("2024-03-20", 48.21, 16.37, 1, 1, 1),
        ("2020-02-16", 78.216667, 15.633333, 1, 1, 1),
        ("2020-04-16", 78.216667, 15.633333, 2, 1, 1),
        ("2020-08-25", 78.216667, 15.633333, 0, 1, 1),
        ("2024-03-18", 90.0, 0.0, 1, 0, 1),
        ("2024-09-20", -90.0, 0.0, 1, 0, 1),
        ("2024-06-21", 90.0, 0.0, 0, 0, 1),
        ("2020-06-10", 0.0, 179.9, 1, 1, 0),
    ] {
        let days = calculator
            .for_date_multiple(
                date.parse().unwrap(),
                &TimeZone::UTC,
                Location {
                    latitude: lat,
                    longitude: lon,
                },
                DELTA_T,
                HORIZONS,
            )
            .unwrap();
        let day = &days[0].1;
        assert_eq!(
            (day.rises.len(), day.sets.len(), day.transits.len()),
            (rises, sets, transits),
            "{date} {lat} {lon}"
        );
        for (horizon, day) in days {
            for (events, direction) in [(&day.rises, 1.0), (&day.sets, -1.0)] {
                for time in events {
                    // Check either side of the 100 ms search bracket.
                    let before = SolarPositions::grena3()
                        .at(
                            time.timestamp() - SignedDuration::from_millis(101),
                            Location {
                                latitude: lat,
                                longitude: lon,
                            },
                            0.0,
                            DELTA_T,
                            None,
                        )
                        .unwrap()
                        .elevation_angle();
                    let after = SolarPositions::grena3()
                        .at(
                            time.timestamp() + SignedDuration::from_millis(101),
                            Location {
                                latitude: lat,
                                longitude: lon,
                            },
                            0.0,
                            DELTA_T,
                            None,
                        )
                        .unwrap()
                        .elevation_angle();
                    assert!(direction * (before - horizon.elevation_angle()) < 0.0);
                    assert!(direction * (after - horizon.elevation_angle()) > 0.0);
                }
            }
            for time in day.transits {
                let time = time.timestamp();
                let before = SolarPositions::grena3()
                    .at(
                        time - SignedDuration::from_millis(101),
                        Location {
                            latitude: 0.0,
                            longitude: lon,
                        },
                        0.0,
                        DELTA_T,
                        None,
                    )
                    .unwrap()
                    .azimuth()
                    .to_radians()
                    .sin();
                let after = SolarPositions::grena3()
                    .at(
                        time + SignedDuration::from_millis(101),
                        Location {
                            latitude: 0.0,
                            longitude: lon,
                        },
                        0.0,
                        DELTA_T,
                        None,
                    )
                    .unwrap()
                    .azimuth()
                    .to_radians()
                    .sin();
                assert!(before > 0.0 && after < 0.0);
            }
        }
    }
}

#[test]
fn finds_seasonal_rise_at_both_poles() {
    for (latitude, start, end) in [
        (90.0, "2024-03-01T00:00:00Z", "2024-04-01T00:00:00Z"),
        (-90.0, "2024-09-01T00:00:00Z", "2024-10-01T00:00:00Z"),
    ] {
        let (start, end) = (instant(start), instant(end));
        let calculator = SolarEvents::new();
        let location = Location {
            latitude,
            longitude: 0.0,
        };
        let horizon = Horizon::SunriseSunset;
        let rise = calculator
            .next_rise(start, end, location, DELTA_T, horizon)
            .unwrap()
            .unwrap();
        for (offset, sign) in [(-101, -1.0), (101, 1.0)] {
            let elevation = SolarPositions::new()
                .at(
                    rise + SignedDuration::from_millis(offset),
                    location,
                    0.0,
                    DELTA_T,
                    None,
                )
                .unwrap()
                .elevation_angle();
            assert!(sign * (elevation - horizon.elevation_angle()) > 0.0);
        }
        assert!(
            calculator
                .next_rise(rise, end, location, DELTA_T, horizon)
                .unwrap()
                .is_none()
        );
        assert!(
            calculator
                .next_set(start, end, location, DELTA_T, horizon)
                .unwrap()
                .is_none()
        );
    }
}

#[test]
fn includes_shallow_endpoint_crossing_without_repeating_it() {
    let root = instant("2024-06-21T12:00:00Z");
    let start = root - SignedDuration::from_secs(10);
    let end = root + SignedDuration::from_secs(10);
    let location = Location {
        latitude: -90.0,
        longitude: 0.0,
    };
    let horizon = Horizon::Custom(
        SolarPositions::grena3()
            .at(root, location, 0.0, 0.0, None)
            .unwrap()
            .elevation_angle(),
    );
    let calculator = SolarEvents::grena3();
    assert_eq!(
        calculator
            .next_rise(start, root, location, 0.0, horizon)
            .unwrap(),
        Some(root)
    );
    assert!(
        calculator
            .next_rise(root, end, location, 0.0, horizon)
            .unwrap()
            .is_none()
    );
}

#[test]
fn assigns_midnight_crossings_to_the_following_date() {
    use solar_positioning::HorizonState;
    let date = date(2024, 3, 20);
    let midnight = TimeZone::UTC.to_timestamp(date.at(0, 0, 0, 0)).unwrap();
    for direction in [-1.0, 1.0] {
        let calculator = SolarEvents::with_provider(
            move |time: JulianDate, _: Location| {
                let phase = std::f64::consts::TAU * (time.julian_date() - 2_460_389.5);
                Ok(EventPosition {
                    elevation: (direction * 0.5 * phase.sin()).asin().to_degrees(),
                    hour_angle: phase.to_degrees(),
                })
            },
            2024..=2024,
        )
        .unwrap();
        let search = |a, b| {
            if direction > 0.0 {
                calculator.next_rise(a, b, ORIGIN, 0.0, Horizon::Custom(0.0))
            } else {
                calculator.next_set(a, b, ORIGIN, 0.0, Horizon::Custom(0.0))
            }
        };
        assert_eq!(
            search(midnight - SignedDuration::from_secs(1), midnight).unwrap(),
            Some(midnight)
        );
        assert!(
            search(midnight, midnight + SignedDuration::from_secs(1))
                .unwrap()
                .is_none()
        );
        for offset in [-1, 0] {
            let day = calculator
                .for_date(
                    if offset < 0 {
                        date.yesterday().unwrap()
                    } else {
                        date
                    },
                    &TimeZone::UTC,
                    ORIGIN,
                    0.0,
                    Horizon::Custom(0.0),
                )
                .unwrap();
            let events = if direction > 0.0 { day.rises } else { day.sets };
            assert_eq!(events.len(), 1);
            assert!(
                (midnight + SignedDuration::from_hours(24 * offset) - events[0].timestamp()).abs()
                    <= SignedDuration::from_millis(100)
            );
            assert_eq!(day.state_at_start, HorizonState::OnHorizon);
        }
    }
}

#[test]
fn uses_continuous_time_across_the_1582_calendar_change() {
    let calculator = SolarEvents::new();
    let mut start = instant("1582-10-04T00:00:00Z");
    let end = instant("1582-10-17T00:00:00Z");
    for day in 4..=16 {
        let transit = calculator
            .next_transit(start, end, 0.0, DELTA_T)
            .unwrap()
            .unwrap();
        assert_eq!(
            TimeZone::UTC.to_datetime(transit).date(),
            date(1582, 10, day)
        );
        start = transit;
    }
    assert!(
        calculator
            .next_transit(start, end, 0.0, DELTA_T)
            .unwrap()
            .is_none()
    );
}

#[test]
fn historical_events_agree_with_positions() {
    let location = Location {
        latitude: 48.21,
        longitude: 16.37,
    };
    let horizon = Horizon::SunriseSunset;
    for date in ["1500-03-20", "1582-10-04", "1582-10-10", "1582-10-15"] {
        let day = SolarEvents::new()
            .for_date(
                date.parse().unwrap(),
                &TimeZone::UTC,
                location,
                DELTA_T,
                horizon,
            )
            .unwrap();
        assert_eq!((day.rises.len(), day.sets.len()), (1, 1));
        for (events, direction) in [(&day.rises, 1.0), (&day.sets, -1.0)] {
            for time in events {
                for offset in [-101, 101] {
                    let position = SolarPositions::new()
                        .at(
                            time.timestamp() + SignedDuration::from_millis(offset),
                            location,
                            0.0,
                            DELTA_T,
                            None,
                        )
                        .unwrap();
                    assert!(
                        direction
                            * (offset as f64)
                            * (position.elevation_angle() - horizon.elevation_angle())
                            > 0.0,
                        "{date}: {time}, offset {offset} ms, elevation {}",
                        position.elevation_angle()
                    );
                }
            }
        }
    }
}

#[test]
fn date_lookback_respects_the_tt_lower_bound() {
    let calculator = SolarEvents::grena3();
    let date = date(2010, 1, 1);
    let zone = offset(-1).to_time_zone();
    let start = zone.to_timestamp(date.at(0, 0, 0, 0)).unwrap();
    let end = start + SignedDuration::from_hours(24);
    for delta_t in [-3600.0, -3599.95, -3599.8] {
        let day = calculator
            .for_date(date, &zone, ORIGIN, delta_t, Horizon::SunriseSunset)
            .unwrap();
        let expected = [
            calculator
                .next_rise(start, end, ORIGIN, delta_t, Horizon::SunriseSunset)
                .unwrap()
                .unwrap(),
            calculator
                .next_transit(start, end, 0.0, delta_t)
                .unwrap()
                .unwrap(),
            calculator
                .next_set(start, end, ORIGIN, delta_t, Horizon::SunriseSunset)
                .unwrap()
                .unwrap(),
        ];
        for (events, time) in [&day.rises, &day.transits, &day.sets]
            .into_iter()
            .zip(expected)
        {
            assert_eq!(events.len(), 1);
            assert!((events[0].timestamp() - time).abs() <= SignedDuration::from_millis(100));
        }
    }
}

#[test]
fn rejects_out_of_range_dates() {
    let calculator = SolarEvents::grena3();
    for date in ["2009-12-31", "2111-01-01"] {
        assert!(
            calculator
                .for_date(
                    date.parse().unwrap(),
                    &TimeZone::UTC,
                    ORIGIN,
                    0.0,
                    Horizon::SunriseSunset
                )
                .is_err()
        );
    }
    let first = instant("2010-01-01T00:00:00Z");
    assert!(
        calculator
            .next_transit(first - SignedDuration::from_nanos(1), first, 0.0, 0.0)
            .is_err()
    );
    let end = instant("2111-01-01T00:00:00Z");
    assert!(
        calculator
            .next_transit(end, end + SignedDuration::from_nanos(1), 0.0, 0.0)
            .is_err()
    );
}

#[test]
fn splitting_searches_preserves_events_for_both_models() {
    let mut seed = 20_260_925_u64;
    let mut random = || {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        (seed >> 32) as u32
    };
    for (calculator, first, years) in [
        (SolarEvents::new(), -1999, 7999),
        (SolarEvents::grena3(), 2011, 99),
    ] {
        for sample in 0..24 {
            let year = first + (random() % years) as i32;
            let start = TimeZone::UTC
                .to_timestamp(date(year as i16, 1, 1).at(0, 0, 0, 0))
                .unwrap()
                + SignedDuration::from_hours(i64::from(random() % 365) * 24)
                + SignedDuration::from_secs(i64::from(random() % 86400));
            let end = start + SignedDuration::from_hours(48);
            let latitude = match sample % 6 {
                0 => 90.0,
                1 => -90.0,
                _ => f64::from(random()) / f64::from(u32::MAX) * 180.0 - 90.0,
            };
            let longitude = f64::from(random()) / f64::from(u32::MAX) * 360.0 - 180.0;
            let location = Location {
                latitude,
                longitude,
            };
            let horizon = if sample % 2 == 0 {
                HORIZONS[(sample / 2) % 4]
            } else {
                Horizon::Custom(f64::from(random()) / f64::from(u32::MAX) * 40.0 - 20.0)
            };
            for kind in 0..3 {
                let search = |a, b| match kind {
                    0 => calculator.next_rise(a, b, location, DELTA_T, horizon),
                    1 => calculator.next_set(a, b, location, DELTA_T, horizon),
                    _ => calculator.next_transit(a, b, longitude, DELTA_T),
                };
                let whole = find_all(start, end, search);
                let mut splits =
                    vec![start + SignedDuration::from_secs(1 + i64::from(random() % 172798))];
                for &event in &whole {
                    for millis in [-101, 0, 101] {
                        splits.push(event + SignedDuration::from_millis(millis));
                    }
                }
                for split in splits {
                    if split <= start || split >= end {
                        continue;
                    }
                    let mut parts = find_all(start, split, search);
                    parts.extend(find_all(split, end, search));
                    compare(&whole, &parts);
                }
                let numeric =
                    find_all_numeric(&calculator, jd(start), jd(end), location, horizon, kind);
                assert_eq!(numeric.len(), whole.len());
                for (&number, &time) in numeric.iter().zip(&whole) {
                    assert!((number - jd(time)).abs() * 86400.0 < 0.1001);
                }
            }
        }
    }
}

fn find_all_numeric(
    calculator: &SolarEvents,
    mut start: f64,
    end: f64,
    location: Location,
    horizon: Horizon,
    kind: u32,
) -> Vec<f64> {
    let mut result = Vec::new();
    loop {
        let next = match kind {
            0 => calculator.next_rise_from_julian(start, end, location, DELTA_T, horizon),
            1 => calculator.next_set_from_julian(start, end, location, DELTA_T, horizon),
            _ => calculator.next_transit_from_julian(start, end, location.longitude, DELTA_T),
        }
        .unwrap();
        let Some(time) = next else {
            return result;
        };
        assert!(time > start && time <= end);
        result.push(time);
        assert!(result.len() < 10);
        start = time;
    }
}
