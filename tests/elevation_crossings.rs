use solar_positioning::{
    ElevationTarget, Error, Horizon, Location, RefractionCorrection, SolarEvents, SolarPositions,
    time::JulianDate,
};

const DELTA_T: f64 = 69.184;
const LOCATION: Location = Location {
    latitude: 52.0,
    longitude: 13.4,
};
const START: f64 = 2_460_389.5;
const END: f64 = START + 1.0;

fn elevation(
    positions: SolarPositions,
    time: f64,
    location: Location,
    target: ElevationTarget,
) -> f64 {
    positions
        .at_from_julian(
            JulianDate::new(time, DELTA_T).unwrap(),
            location,
            target.height,
            target.refraction,
        )
        .unwrap()
        .elevation_angle()
}

fn crossing(
    positions: SolarPositions,
    rising: bool,
    start: f64,
    end: f64,
    location: Location,
    target: ElevationTarget,
) -> Option<f64> {
    if rising {
        positions.next_elevation_rise_from_julian(start, end, location, DELTA_T, target)
    } else {
        positions.next_elevation_set_from_julian(start, end, location, DELTA_T, target)
    }
    .unwrap()
}

fn standard_crossing(events: SolarEvents, rising: bool, start: f64, end: f64, angle: f64) -> f64 {
    if rising {
        events.next_rise_from_julian(start, end, LOCATION, DELTA_T, Horizon::Custom(angle))
    } else {
        events.next_set_from_julian(start, end, LOCATION, DELTA_T, Horizon::Custom(angle))
    }
    .unwrap()
    .unwrap()
}

fn crosses(before: f64, after: f64, angle: f64, rising: bool) -> bool {
    if rising {
        before < angle && after > angle
    } else {
        before > angle && after < angle
    }
}

fn assert_transitions(
    positions: SolarPositions,
    start: f64,
    location: Location,
    target: ElevationTarget,
) {
    for rising in [true, false] {
        if let Some(event) = crossing(positions, rising, start, start + 1.0, location, target) {
            assert!(event > start && event <= start + 1.0);
            let before = elevation(positions, event - 0.1 / 86400.0, location, target);
            let after = elevation(positions, event + 0.1 / 86400.0, location, target);
            assert!(
                crosses(before, after, target.angle, rising),
                "{positions:?}, {location:?}, {target:?}: {before} -> {after}"
            );
        } else {
            let mut before = elevation(positions, start, location, target);
            for i in 1..=96 {
                let after = elevation(positions, start + f64::from(i) / 96.0, location, target);
                assert!(
                    !crosses(before, after, target.angle, rising),
                    "missing crossing: {positions:?}, {location:?}, {target:?}"
                );
                before = after;
            }
        }
    }
}

#[test]
fn crossings_match_position_predicates_across_settings_and_seasons() {
    for (positions, heights) in [
        (SolarPositions::new(), &[0.0, 2500.0][..]),
        (SolarPositions::grena3(), &[0.0][..]),
    ] {
        for month in [1, 3, 6] {
            let start = JulianDate::from_utc(2024, month, 21, 0, 0, 0.0, DELTA_T)
                .unwrap()
                .julian_date();
            for latitude in [-35.0, 0.0, 52.0] {
                for &height in heights {
                    for refraction in [
                        None,
                        Some(RefractionCorrection::standard()),
                        Some(RefractionCorrection::new(850.0, -10.0).unwrap()),
                        Some(RefractionCorrection::new(1100.0, -70.0).unwrap()),
                    ] {
                        for angle in [-10.0, 0.0, 10.0, 40.0] {
                            assert_transitions(
                                positions,
                                start,
                                Location {
                                    latitude,
                                    ..LOCATION
                                },
                                ElevationTarget {
                                    angle,
                                    height,
                                    refraction,
                                },
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn spa_crossings_use_observer_height() {
    let positions = SolarPositions::new();
    let location = Location {
        latitude: 63.0,
        ..LOCATION
    };
    let target = ElevationTarget {
        angle: 0.0,
        height: 5000.0,
        refraction: None,
    };
    let start = JulianDate::from_utc(2024, 12, 21, 0, 0, 0.0, DELTA_T)
        .unwrap()
        .julian_date();
    // Shallow winter crossings make the height effect exceed the 1 ms search resolution.
    for rising in [true, false] {
        let event = crossing(positions, rising, start, start + 1.0, location, target).unwrap();
        let before = elevation(positions, event - 0.0011 / 86400.0, location, target);
        let after = elevation(positions, event + 0.0011 / 86400.0, location, target);
        assert!(
            crosses(before, after, target.angle, rising),
            "height-adjusted crossing: {before} -> {after}, rising={rising}"
        );
    }
}

#[test]
fn unrefracted_sea_level_crossings_match_existing_events_and_resume() {
    for (positions, events) in [
        (SolarPositions::new(), SolarEvents::new()),
        (SolarPositions::grena3(), SolarEvents::grena3()),
    ] {
        for angle in [-18.0, 0.0, 10.0] {
            let target = ElevationTarget {
                angle,
                height: 0.0,
                refraction: None,
            };
            for rising in [true, false] {
                let event = crossing(positions, rising, START, END, LOCATION, target).unwrap();
                assert!(
                    (event - standard_crossing(events, rising, START, END, angle)).abs() * 86400.0
                        < 0.0011
                );
                assert!(crossing(positions, rising, event, END, LOCATION, target).is_none());
                let bounded = crossing(positions, rising, START, event, LOCATION, target).unwrap();
                assert!((bounded - event).abs() * 86400.0 < 0.0011);
            }
        }
    }
}

#[test]
fn thresholds_inside_refraction_jumps_cross_at_the_geometric_cutoff() {
    for (positions, events, cutoff, angle) in [
        (SolarPositions::new(), SolarEvents::new(), -0.83337, -0.5),
        (SolarPositions::grena3(), SolarEvents::grena3(), 0.0, 0.2),
    ] {
        let target = ElevationTarget {
            angle,
            height: 0.0,
            refraction: Some(RefractionCorrection::standard()),
        };
        for rising in [true, false] {
            let event = crossing(positions, rising, START, END, LOCATION, target).unwrap();
            assert!(
                (event - standard_crossing(events, rising, START, END, cutoff)).abs() * 86400.0
                    < 0.0011
            );
        }
        assert_transitions(positions, START, LOCATION, target);
    }
}

#[test]
fn missing_crossings_and_invalid_inputs_remain_distinct() {
    let target = ElevationTarget {
        angle: 0.0,
        height: 0.0,
        refraction: None,
    };
    for positions in [SolarPositions::new(), SolarPositions::grena3()] {
        let query = |target| {
            positions.next_elevation_rise_from_julian(START, END, LOCATION, DELTA_T, target)
        };
        let polar = Location {
            latitude: 90.0,
            longitude: 0.0,
        };
        for rising in [true, false] {
            assert!(crossing(positions, rising, 2_460_482.5, 2_460_483.5, polar, target).is_none());
            assert!(crossing(positions, rising, START, START, LOCATION, target).is_none());
        }
        for angle in [f64::NAN, -91.0, 91.0] {
            assert!(matches!(
                query(ElevationTarget { angle, ..target }),
                Err(Error::InvalidElevationAngle { .. })
            ));
        }
        for height in [f64::NAN, f64::INFINITY] {
            assert!(matches!(
                query(ElevationTarget { height, ..target }),
                Err(Error::InvalidHeight { .. })
            ));
        }
        for (a, b, delta_t) in [
            (END, START, DELTA_T),
            (START, END, f64::NAN),
            (f64::NAN, END, DELTA_T),
            (0.0, 1.0, DELTA_T),
        ] {
            assert!(matches!(
                positions.next_elevation_rise_from_julian(a, b, LOCATION, delta_t, target),
                Err(Error::InvalidDateTime { .. })
            ));
        }
        assert!(matches!(
            positions.next_elevation_rise_from_julian(
                START,
                END,
                Location {
                    latitude: 91.0,
                    ..LOCATION
                },
                DELTA_T,
                target
            ),
            Err(Error::InvalidLatitude { .. })
        ));
        let unreachable = ElevationTarget {
            angle: 90.0,
            refraction: Some(RefractionCorrection::standard()),
            ..target
        };
        assert!(query(unreachable).unwrap().is_none());
        assert!(
            positions
                .next_elevation_rise_from_julian(END, START, LOCATION, DELTA_T, unreachable)
                .is_err()
        );
        let excessive = ElevationTarget {
            refraction: Some(RefractionCorrection::new(2000.0, -250.0).unwrap()),
            ..target
        };
        assert!(matches!(
            query(excessive),
            Err(Error::ComputationError { .. })
        ));
    }
    assert!(matches!(
        SolarPositions::grena3().next_elevation_rise_from_julian(
            START,
            END,
            LOCATION,
            DELTA_T,
            ElevationTarget {
                height: 1.0,
                ..target
            }
        ),
        Err(Error::InvalidHeight { .. })
    ));
}

#[cfg(feature = "chrono")]
#[test]
fn chrono_matches_numeric_search_and_preserves_dst_offsets() {
    use chrono::{Days, Duration, TimeZone};
    use chrono_tz::Europe::Berlin;
    for (month, day, offset) in [(3, 31, 7200), (10, 27, 3600)] {
        let start = Berlin.with_ymd_and_hms(2024, month, day, 0, 0, 0).unwrap();
        let end = start.checked_add_days(Days::new(1)).unwrap();
        let a = JulianDate::from_datetime(&start, DELTA_T)
            .unwrap()
            .julian_date();
        let b = JulianDate::from_datetime(&end, DELTA_T)
            .unwrap()
            .julian_date();
        for (positions, height) in [
            (SolarPositions::new(), 1500.0),
            (SolarPositions::grena3(), 0.0),
        ] {
            let target = ElevationTarget {
                angle: 10.0,
                height,
                refraction: Some(RefractionCorrection::standard()),
            };
            for (rising, result) in [
                (
                    true,
                    positions.next_elevation_rise(&start, &end, LOCATION, DELTA_T, target),
                ),
                (
                    false,
                    positions.next_elevation_set(&start, &end, LOCATION, DELTA_T, target),
                ),
            ] {
                let event = result.unwrap().unwrap();
                let numeric = crossing(positions, rising, a, b, LOCATION, target).unwrap();
                assert!(
                    (JulianDate::from_datetime(&event, DELTA_T)
                        .unwrap()
                        .julian_date()
                        - numeric)
                        .abs()
                        * 86400.0
                        < 0.0011
                );
                assert_eq!(event.fixed_offset().offset().local_minus_utc(), offset);
                let before = positions
                    .at(
                        &(event - Duration::milliseconds(100)),
                        LOCATION,
                        height,
                        DELTA_T,
                        target.refraction,
                    )
                    .unwrap()
                    .elevation_angle();
                let after = positions
                    .at(
                        &(event + Duration::milliseconds(100)),
                        LOCATION,
                        height,
                        DELTA_T,
                        target.refraction,
                    )
                    .unwrap()
                    .elevation_angle();
                assert!(crosses(before, after, target.angle, rising));
            }
        }
    }
}
