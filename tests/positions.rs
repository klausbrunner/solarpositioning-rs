use solar_positioning::{
    Error, Location, PreparedPositions, RefractionCorrection, SolarPositions, time::JulianDate,
};

const LOCATION: Location = Location {
    latitude: 48.21,
    longitude: 16.37,
};

#[test]
fn prepared_positions_match_single_calls_and_remain_independent() {
    fn assert_traits<T: Copy + Send + Sync>() {}
    assert_traits::<SolarPositions>();
    assert_traits::<PreparedPositions>();
    for (positions, height) in [
        (SolarPositions::default(), 190.0),
        (SolarPositions::grena3(), 0.0),
    ] {
        for month in [1, 3, 6, 9, 12] {
            let time = JulianDate::from_utc(2024, month, 21, 12, 30, 12.5, 69.184).unwrap();
            let prepared = positions.for_time_from_julian(time);
            let original = prepared.at(LOCATION, height, None).unwrap();
            let later = JulianDate::new(time.julian_date() + 0.25, time.delta_t() + 1.0).unwrap();
            let other = positions.for_time_from_julian(later);
            assert_eq!(
                other.at(LOCATION, height, None),
                positions.at_from_julian(later, LOCATION, height, None)
            );
            assert_eq!(original, prepared.at(LOCATION, height, None).unwrap());
            for latitude in [-90.0, -48.21, 0.0, 48.21, 90.0] {
                for longitude in [-180.0, 0.0, 16.37, 180.0] {
                    let location = Location {
                        latitude,
                        longitude,
                    };
                    for refraction in [None, Some(RefractionCorrection::standard())] {
                        let expected = positions
                            .at_from_julian(time, location, height, refraction)
                            .unwrap();
                        assert_eq!(expected, prepared.at(location, height, refraction).unwrap());
                    }
                }
            }
        }
    }
}

#[test]
fn published_position_examples() {
    let time = JulianDate::from_utc(2003, 10, 17, 19, 30, 30.0, 67.0).unwrap();
    let position = SolarPositions::new()
        .at_from_julian(
            time,
            Location {
                latitude: 39.742476,
                longitude: -105.1786,
            },
            1830.14,
            Some(RefractionCorrection::new(820.0, 11.0).unwrap()),
        )
        .unwrap();
    assert!((position.azimuth() - 194.340241).abs() < 0.000001);
    assert!((position.zenith_angle() - 50.111622).abs() < 0.000001);

    // Grena's sample reports angles to five decimal places in radians.
    let time = JulianDate::from_utc(2012, 1, 1, 11, 0, 0.0, 65.0).unwrap();
    let position = SolarPositions::grena3()
        .at_from_julian(
            time,
            Location {
                latitude: 0.73117_f64.to_degrees(),
                longitude: 0.21787_f64.to_degrees(),
            },
            0.0,
            Some(RefractionCorrection::new(1000.0, 20.0).unwrap()),
        )
        .unwrap();
    assert!(
        (position.azimuth() - (-0.0591845 + core::f64::consts::PI).to_degrees()).abs() < 0.0004
    );
    assert!((position.zenith_angle() - 1.13381_f64.to_degrees()).abs() < 0.0004);
}

#[test]
fn spa_handles_zenith_and_nadir_roundoff() {
    // These zenith/nadir cases round the cosine just outside [-1, 1].
    for (jd, latitude, longitude, elevation) in [
        (
            2_460_316.0,
            -22.521_369_002_488_342,
            1.401_502_537_686_724_3,
            90.0,
        ),
        (
            2_460_320.0,
            21.989_156_589_016_382,
            -178.172_550_684_033_8,
            -90.0,
        ),
    ] {
        let position = SolarPositions::new()
            .at_from_julian(
                JulianDate::new(jd, 69.184).unwrap(),
                Location {
                    latitude,
                    longitude,
                },
                0.0,
                None,
            )
            .unwrap();
        assert!((position.elevation_angle() - elevation).abs() < 0.000_02);
        assert!(position.azimuth().is_finite());
    }
}

#[test]
fn rejects_invalid_coordinates_and_heights() {
    let time = JulianDate::from_utc(2024, 6, 21, 12, 0, 0.0, 69.0).unwrap();
    for positions in [SolarPositions::new(), SolarPositions::grena3()] {
        let prepared = positions.for_time_from_julian(time);
        for location in [
            Location {
                latitude: 91.0,
                longitude: 0.0,
            },
            Location {
                latitude: 0.0,
                longitude: -181.0,
            },
            Location {
                latitude: f64::NAN,
                longitude: 0.0,
            },
            Location {
                latitude: 0.0,
                longitude: f64::INFINITY,
            },
        ] {
            assert!(prepared.at(location, 0.0, None).is_err());
            assert!(positions.at_from_julian(time, location, 0.0, None).is_err());
        }
        for height in [f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
            assert!(matches!(
                prepared.at(LOCATION, height, None),
                Err(Error::InvalidHeight { .. })
            ));
        }
    }
    for height in [-50.0, 190.0] {
        assert!(
            SolarPositions::new()
                .at_from_julian(time, LOCATION, height, None)
                .is_ok()
        );
        assert_eq!(
            SolarPositions::grena3().at_from_julian(time, LOCATION, height, None),
            Err(Error::InvalidHeight { value: height })
        );
    }
}

#[cfg(feature = "chrono")]
#[test]
fn chrono_matches_numeric_time_across_offsets_and_year_boundaries() {
    use chrono::{DateTime, Duration, FixedOffset, Utc};
    for text in [
        "2024-06-21T08:30:12.345678+02:00",
        "2024-01-01T00:00:00.123456789+14:00",
        "2024-12-31T23:59:59.987654321-11:00",
    ] {
        let time = text.parse::<DateTime<FixedOffset>>().unwrap();
        let jd = JulianDate::from_datetime(&time, 69.184).unwrap();
        for (positions, height) in [
            (SolarPositions::new(), -50.0),
            (SolarPositions::grena3(), 0.0),
        ] {
            for refraction in [None, Some(RefractionCorrection::standard())] {
                let expected = positions
                    .at_from_julian(jd, LOCATION, height, refraction)
                    .unwrap();
                assert_eq!(
                    expected,
                    positions
                        .at(&time, LOCATION, height, 69.184, refraction)
                        .unwrap()
                );
                assert_eq!(
                    expected,
                    positions
                        .at(
                            &time.with_timezone(&Utc),
                            LOCATION,
                            height,
                            69.184,
                            refraction
                        )
                        .unwrap()
                );
                assert_eq!(
                    expected,
                    positions
                        .for_time(&time, 69.184)
                        .unwrap()
                        .at(LOCATION, height, refraction)
                        .unwrap()
                );
                assert_ne!(
                    expected,
                    positions
                        .at(
                            &(time + Duration::milliseconds(500)),
                            LOCATION,
                            height,
                            69.184,
                            refraction
                        )
                        .unwrap()
                );
            }
            for delta_t in [f64::NAN, f64::INFINITY] {
                assert!(positions.for_time(&time, delta_t).is_err());
                assert!(
                    positions
                        .at(&time, LOCATION, height, delta_t, None)
                        .is_err()
                );
            }
        }
    }
}
