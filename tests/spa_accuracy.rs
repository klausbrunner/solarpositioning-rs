//! SPA accuracy tests using reference data from NREL's reference implementation.

use chrono::{DateTime, FixedOffset, TimeZone, Utc};
use csv::ReaderBuilder;
use solar_positioning::{RefractionCorrection, spa};
use std::error::Error;
use std::fs::File;

#[allow(clippy::type_complexity)]
fn load_spa_reference_data() -> Result<Vec<(DateTime<Utc>, f64, f64, f64, f64)>, Box<dyn Error>> {
    let file = File::open("tests/data/spa_reference_testdata.csv")?;
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(file);

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result?;
        if record.len() >= 5 {
            let datetime = record[0].parse::<DateTime<Utc>>()?;
            let latitude: f64 = record[1].parse()?;
            let longitude: f64 = record[2].parse()?;
            let expected_azimuth: f64 = record[3].parse()?;
            let expected_zenith: f64 = record[4].parse()?;
            records.push((
                datetime,
                latitude,
                longitude,
                expected_azimuth,
                expected_zenith,
            ));
        }
    }

    Ok(records)
}

#[test]
fn test_spa_accuracy_against_nrel_reference() -> Result<(), Box<dyn Error>> {
    let test_records = load_spa_reference_data()?;

    assert!(!test_records.is_empty(), "No test records loaded");

    println!("Testing {} SPA reference records", test_records.len());

    let mut max_azimuth_error = 0.0_f64;
    let mut max_zenith_error = 0.0_f64;
    let mut error_count = 0;

    // Test parameters from NREL reference: elevation=0, pressure=1000, temperature=10, deltaT=0
    for (i, (datetime, latitude, longitude, expected_azimuth, expected_zenith)) in
        test_records.iter().enumerate()
    {
        let datetime_fixed = FixedOffset::east_opt(0)
            .unwrap()
            .from_utc_datetime(&datetime.naive_utc());
        match spa::solar_position(
            datetime_fixed,
            *latitude,
            *longitude,
            0.0,                                                    // elevation (meters)
            0.0, // deltaT (seconds) - reference uses 0
            Some(RefractionCorrection::new(1000.0, 10.0).unwrap()), // atmospheric conditions
        ) {
            Ok(position) => {
                let azimuth_error = (position.azimuth() - expected_azimuth).abs();
                let zenith_error = (position.zenith_angle() - expected_zenith).abs();

                max_azimuth_error = max_azimuth_error.max(azimuth_error);
                max_zenith_error = max_zenith_error.max(zenith_error);

                // SPA should be accurate to ±0.0003 degrees according to NREL
                // We'll use a more relaxed tolerance initially to debug any issues
                let tolerance = 0.001; // 0.001 degrees = 3.6 arcseconds

                if azimuth_error > tolerance || zenith_error > tolerance {
                    println!(
                        "Record {}: DateTime={}, Lat={:.6}, Lon={:.6}",
                        i, datetime, latitude, longitude
                    );
                    println!(
                        "  Expected: Az={:.6}°, Zen={:.6}°",
                        expected_azimuth, expected_zenith
                    );
                    println!(
                        "  Actual:   Az={:.6}°, Zen={:.6}°",
                        position.azimuth(),
                        position.zenith_angle()
                    );
                    println!(
                        "  Error:    Az={:.6}°, Zen={:.6}°",
                        azimuth_error, zenith_error
                    );

                    error_count += 1;
                    if error_count > 10 {
                        println!("... (showing first 10 errors only)");
                        break;
                    }
                }
            }
            Err(e) => {
                panic!("SPA calculation failed for record {}: {}", i, e);
            }
        }
    }

    println!("Maximum azimuth error: {:.6}°", max_azimuth_error);
    println!("Maximum zenith error: {:.6}°", max_zenith_error);

    // SPA algorithm should achieve NREL specification of ±0.0003 degrees
    let nrel_tolerance = 0.0003; // 0.0003 degrees = 1.08 arcseconds (NREL specification)

    assert!(
        max_azimuth_error < nrel_tolerance,
        "Maximum azimuth error {:.6}° exceeds tolerance {:.6}°",
        max_azimuth_error,
        nrel_tolerance
    );

    assert!(
        max_zenith_error < nrel_tolerance,
        "Maximum zenith error {:.6}° exceeds tolerance {:.6}°",
        max_zenith_error,
        nrel_tolerance
    );

    Ok(())
}

#[test]
fn test_spa_specific_cases() {
    // Test a few specific cases that are known to be correct

    // Test case: 2000-01-01 12:00 UTC at Greenwich (lat=0, lon=0)
    let utc_datetime = "2000-01-01T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let datetime = FixedOffset::east_opt(0)
        .unwrap()
        .from_utc_datetime(&utc_datetime.naive_utc());
    let result = spa::solar_position(
        datetime,
        0.0,                                    // Greenwich latitude
        0.0,                                    // Greenwich longitude
        0.0,                                    // sea level
        63.86,                                  // deltaT for year 2000 (from DeltaT estimation)
        Some(RefractionCorrection::standard()), // standard atmospheric conditions
    );

    assert!(result.is_ok());
    let position = result.unwrap();

    // At Greenwich on Jan 1, 2000 at solar noon, sun should be roughly south
    // (exact values depend on equation of time, but azimuth should be close to 180°)
    assert!(
        position.azimuth() > 170.0 && position.azimuth() < 190.0,
        "Azimuth {:.3}° not near south at Greenwich solar noon",
        position.azimuth()
    );

    // Zenith angle should be approximately equal to latitude + solar declination
    // On Jan 1, solar declination is about -23°, so zenith should be around 23°
    assert!(
        position.zenith_angle() > 20.0 && position.zenith_angle() < 26.0,
        "Zenith angle {:.3}° not reasonable for Greenwich on Jan 1",
        position.zenith_angle()
    );
}
