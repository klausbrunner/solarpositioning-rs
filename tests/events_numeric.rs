//! Core search guarantees, including builds without chrono or allocation.
use solar_positioning::SolarPositions;
use solar_positioning::{Error, EventPosition, Horizon, Location, SolarEvents, time::JulianDate};

const ORIGIN: Location = Location {
    latitude: 0.0,
    longitude: 0.0,
};

fn find_all(
    mut start: f64,
    end: f64,
    next: impl Fn(f64, f64) -> solar_positioning::Result<Option<f64>>,
) -> Vec<f64> {
    let mut results = Vec::new();
    while let Some(time) = next(start, end).unwrap() {
        assert!(time > start && time <= end, "{start} < {time} <= {end}");
        results.push(time);
        assert!(results.len() < 10, "repeated crossing");
        start = time;
    }
    results
}

#[test]
fn finds_known_crossings_and_rejects_tangencies() {
    // sin(elevation) = -0.5*sin(phase/2)^2 touches zero at noon.
    // Negative horizons have two analytically known crossings, including a ~1 s pair.
    let calculator = SolarEvents::with_provider(
        |time: JulianDate, _: Location| {
            let phase = std::f64::consts::TAU * (time.julian_date() - 2_460_390.0);
            Ok(EventPosition {
                elevation: (-0.5 * (phase / 2.0).sin().powi(2)).asin().to_degrees(),
                hour_angle: phase.to_degrees(),
            })
        },
        2024..=2024,
    )
    .unwrap();
    for horizon in [-10.0_f64, -1e-8, 0.0, 1e-8] {
        let rises = find_all(2_460_389.5, 2_460_390.5, |a, b| {
            calculator.next_rise_from_julian(a, b, ORIGIN, 0.0, Horizon::Custom(horizon))
        });
        let sets = find_all(2_460_389.5, 2_460_390.5, |a, b| {
            calculator.next_set_from_julian(a, b, ORIGIN, 0.0, Horizon::Custom(horizon))
        });
        if horizon >= 0.0 {
            assert!(rises.is_empty() && sets.is_empty());
        } else {
            let days = (-2.0 * horizon.to_radians().sin()).sqrt().asin() / std::f64::consts::PI;
            assert_eq!(rises.len(), 1);
            assert_eq!(sets.len(), 1);
            assert!((rises[0] - (2_460_390.0 - days)).abs() * 86400.0 < 0.0011);
            assert!((sets[0] - (2_460_390.0 + days)).abs() * 86400.0 < 0.0011);
        }
    }
}

#[test]
fn preserves_bounds_and_does_not_repeat_returned_events() {
    let calculator = SolarEvents::new();
    let (start, end) = (2_460_389.5, 2_460_390.5);
    let transit = calculator
        .next_transit_from_julian(start, end, 0.0, 69.184)
        .unwrap()
        .unwrap();
    let bounded = calculator
        .next_transit_from_julian(start, transit, 0.0, 69.184)
        .unwrap()
        .unwrap();
    assert!((bounded - transit).abs() * 86400.0 < 0.0011);
    assert!(
        calculator
            .next_transit_from_julian(transit, end, 0.0, 69.184)
            .unwrap()
            .is_none()
    );
    assert!(
        calculator
            .next_transit_from_julian(start, transit - 0.002 / 86400.0, 0.0, 69.184)
            .unwrap()
            .is_none()
    );
    assert!(
        calculator
            .next_transit_from_julian(start, start, 0.0, 69.184)
            .unwrap()
            .is_none()
    );
    for longitude in [-179.9, 0.0, 179.9] {
        let location = Location {
            latitude: 0.0,
            longitude,
        };
        assert_eq!(
            find_all(start, end, |a, b| calculator.next_rise_from_julian(
                a,
                b,
                location,
                69.184,
                Horizon::SunriseSunset
            ))
            .len(),
            1
        );
        assert_eq!(
            find_all(start, end, |a, b| calculator.next_set_from_julian(
                a,
                b,
                location,
                69.184,
                Horizon::SunriseSunset
            ))
            .len(),
            1
        );
    }
}

#[test]
fn custom_provider_receives_time_coordinates_and_delta_t() {
    let calculator = SolarEvents::with_provider(
        |time: JulianDate, location: Location| {
            let hour_angle =
                (time.julian_ephemeris_day() - 2_460_390.0) * 360.0 + location.longitude;
            Ok(EventPosition {
                elevation: (location.latitude.to_radians().cos() * hour_angle.to_radians().cos())
                    .asin()
                    .to_degrees(),
                hour_angle,
            })
        },
        2024..=2024,
    )
    .unwrap();
    let location = Location {
        latitude: 30.0,
        longitude: 15.0,
    };
    let (start, end) = (2_460_389.5, 2_460_390.5);
    for (time, hour) in [
        (
            calculator
                .next_rise_from_julian(start, end, location, 60.0, Horizon::Custom(0.0))
                .unwrap()
                .unwrap(),
            5.0,
        ),
        (
            calculator
                .next_transit_from_julian(start, end, 15.0, 60.0)
                .unwrap()
                .unwrap(),
            11.0,
        ),
        (
            calculator
                .next_set_from_julian(start, end, location, 60.0, Horizon::Custom(0.0))
                .unwrap()
                .unwrap(),
            17.0,
        ),
    ] {
        assert!(((time - start) * 86400.0 - (hour * 3600.0 - 60.0)).abs() < 0.0011);
    }
    assert!(
        calculator
            .next_transit_from_julian(start - 366.0, end - 366.0, 0.0, 0.0)
            .is_err()
    );
}

#[test]
fn rejects_invalid_inputs_and_provider_results() {
    let calculator = SolarEvents::new();
    let (start, end) = (2_460_389.5, 2_460_390.5);
    for (a, b, dt) in [
        (end, start, 0.0),
        (start, end, f64::NAN),
        (f64::NAN, end, 0.0),
        (start, f64::INFINITY, 0.0),
        (1.0, end, 0.0),
    ] {
        assert!(calculator.next_transit_from_julian(a, b, 0.0, dt).is_err());
    }
    for (lat, lon, horizon) in [
        (91.0, 0.0, 0.0),
        (0.0, 181.0, 0.0),
        (0.0, 0.0, 91.0),
        (f64::NAN, 0.0, 0.0),
        (0.0, 0.0, f64::NAN),
    ] {
        assert!(
            calculator
                .next_rise_from_julian(
                    start,
                    end,
                    Location {
                        latitude: lat,
                        longitude: lon
                    },
                    0.0,
                    Horizon::Custom(horizon)
                )
                .is_err()
        );
    }
    let fail = Error::ComputationError {
        message: "provider failed",
    };
    let provider = |_: JulianDate, _: Location| Err(fail.clone());
    let custom = SolarEvents::with_provider(provider, 2024..=2024).unwrap();
    assert_eq!(
        custom.next_transit_from_julian(start, end, 0.0, 0.0),
        Err(fail.clone())
    );
    for (first, last) in [(-2001, 2024), (2024, 6001), (2024, 2023)] {
        assert!(SolarEvents::with_provider(provider, first..=last).is_err());
    }
    for (elevation, hour_angle) in [(f64::NAN, 0.0), (0.0, f64::INFINITY)] {
        let custom = SolarEvents::with_provider(
            move |_, _| {
                Ok(EventPosition {
                    elevation,
                    hour_angle,
                })
            },
            2024..=2024,
        )
        .unwrap();
        assert!(
            custom
                .next_transit_from_julian(start, end, 0.0, 0.0)
                .is_err()
        );
    }
    assert!(JulianDate::new(f64::NAN, 0.0).is_err());
    assert!(JulianDate::new(start, f64::INFINITY).is_err());
}

#[test]
fn grena_search_handles_zenith_roundoff() {
    let start = JulianDate::from_utc(2024, 1, 13, 0, 0, 0.0, 69.184)
        .unwrap()
        .julian_date();
    let end = start + 1.0;
    // These coordinates put the Sun at the zenith at the search midpoint.
    let location = Location {
        latitude: -21.51229285091162,
        longitude: 2.1219181064613295,
    };
    let horizon = Horizon::SunriseSunset;
    let rise = SolarEvents::grena3()
        .next_rise_from_julian(start, end, location, 69.184, horizon)
        .unwrap()
        .unwrap();
    assert!(rise > start && rise < end);
    for direction in [-1.0, 1.0] {
        let elevation = SolarPositions::grena3()
            .at_from_julian(
                JulianDate::new(rise + direction * 0.002 / 86400.0, 69.184).unwrap(),
                location,
                0.0,
                None,
            )
            .unwrap()
            .elevation_angle();
        assert!(direction * (elevation - horizon.elevation_angle()) > 0.0);
    }
}

#[test]
fn provider_ranges_apply_to_both_ut_and_tt() {
    let calculator = SolarEvents::grena3();
    let start = JulianDate::from_utc(2010, 1, 1, 0, 0, 0.0, 0.0)
        .unwrap()
        .julian_date();
    let end = JulianDate::from_utc(2111, 1, 1, 0, 0, 0.0, 0.0)
        .unwrap()
        .julian_date();
    assert!(
        calculator
            .next_transit_from_julian(start, start + 1.0, 0.0, 0.0)
            .unwrap()
            .is_some()
    );
    assert!(
        calculator
            .next_transit_from_julian(end - 1.0, end, 0.0, 0.0)
            .unwrap()
            .is_some()
    );
    for (a, b, dt) in [
        (start - 1.0, start, 0.0),
        (end, end + 1.0, 0.0),
        (start, start + 1.0, -1.0),
        (end - 1.0, end, 1.0),
    ] {
        assert!(calculator.next_transit_from_julian(a, b, 0.0, dt).is_err());
    }
}

#[test]
#[cfg(feature = "alloc")]
fn numeric_collection_uses_half_open_intervals_and_independent_states() {
    use solar_positioning::HorizonState;
    let calculator = SolarEvents::new();
    let location = Location {
        latitude: 70.0,
        longitude: 0.0,
    };
    let start = JulianDate::from_utc(2024, 12, 21, 0, 0, 0.0, 69.184)
        .unwrap()
        .julian_date();
    let night = calculator
        .for_interval_from_julian(start, start + 1.0, location, 69.184, Horizon::SunriseSunset)
        .unwrap();
    let twilight = calculator
        .for_interval_from_julian(start, start + 1.0, location, 69.184, Horizon::CivilTwilight)
        .unwrap();
    assert!(night.always_below());
    assert!(!twilight.always_below());
    assert_eq!(twilight.state_at_start, HorizonState::Below);
    assert_eq!(twilight.rises.len(), 1);
    assert_eq!(night.transits, twilight.transits);
    for horizon in [Horizon::Custom(-90.0), Horizon::Custom(90.0)] {
        let day = calculator
            .for_interval_from_julian(start, start + 1.0, ORIGIN, 69.184, horizon)
            .unwrap();
        assert_eq!(day.always_above(), horizon.elevation_angle() < 0.0);
        assert_eq!(day.always_below(), horizon.elevation_angle() > 0.0);
    }
    let empty = calculator
        .for_interval_from_julian(start, start, ORIGIN, 69.184, Horizon::SunriseSunset)
        .unwrap();
    assert!(empty.rises.is_empty() && empty.sets.is_empty() && empty.transits.is_empty());
    assert!(!empty.always_above() && !empty.always_below());
}
