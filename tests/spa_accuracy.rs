//! SPA accuracy tests using reference data from NREL's reference implementation.

use chrono::{DateTime, Utc};
use csv::ReaderBuilder;
use solar_positioning::spa;
use std::error::Error;
use std::fs::File;

#[derive(Debug)]
struct SpaTestRecord {
    datetime: DateTime<Utc>,
    latitude: f64,
    longitude: f64,
    expected_azimuth: f64,
    expected_zenith: f64,
}

impl SpaTestRecord {
    fn from_csv_record(record: &csv::StringRecord) -> Result<Self, Box<dyn Error>> {
        let datetime_str = record.get(0).ok_or("Missing datetime")?;
        let latitude: f64 = record.get(1).ok_or("Missing latitude")?.parse()?;
        let longitude: f64 = record.get(2).ok_or("Missing longitude")?.parse()?;
        let expected_azimuth: f64 = record.get(3).ok_or("Missing azimuth")?.parse()?;
        let expected_zenith: f64 = record.get(4).ok_or("Missing zenith")?.parse()?;

        let datetime = datetime_str.parse::<DateTime<Utc>>()?;

        Ok(Self {
            datetime,
            latitude,
            longitude,
            expected_azimuth,
            expected_zenith,
        })
    }
}

fn load_spa_reference_data() -> Result<Vec<SpaTestRecord>, Box<dyn Error>> {
    let file = File::open("tests/data/spa_reference_testdata.csv")?;
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(file);

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result?;
        if !record.is_empty() && record.len() >= 5 {
            match SpaTestRecord::from_csv_record(&record) {
                Ok(test_record) => records.push(test_record),
                Err(e) => {
                    eprintln!("Warning: Failed to parse record {:?}: {}", record, e);
                    continue;
                }
            }
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
    for (i, record) in test_records.iter().enumerate() {
        match spa::solar_position(
            record.datetime,
            record.latitude,
            record.longitude,
            0.0,    // elevation (meters)
            0.0,    // deltaT (seconds) - reference uses 0
            1000.0, // pressure (millibars)
            10.0,   // temperature (°C)
        ) {
            Ok(position) => {
                let azimuth_error = (position.azimuth() - record.expected_azimuth).abs();
                let zenith_error = (position.zenith_angle() - record.expected_zenith).abs();

                max_azimuth_error = max_azimuth_error.max(azimuth_error);
                max_zenith_error = max_zenith_error.max(zenith_error);

                // SPA should be accurate to ±0.0003 degrees according to NREL
                // We'll use a more relaxed tolerance initially to debug any issues
                let tolerance = 0.001; // 0.001 degrees = 3.6 arcseconds

                if azimuth_error > tolerance || zenith_error > tolerance {
                    println!(
                        "Record {}: DateTime={}, Lat={:.6}, Lon={:.6}",
                        i, record.datetime, record.latitude, record.longitude
                    );
                    println!(
                        "  Expected: Az={:.6}°, Zen={:.6}°",
                        record.expected_azimuth, record.expected_zenith
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
    let datetime = "2000-01-01T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let result = spa::solar_position(
        datetime, 0.0,     // Greenwich latitude
        0.0,     // Greenwich longitude
        0.0,     // sea level
        63.86,   // deltaT for year 2000 (from DeltaT estimation)
        1013.25, // standard pressure
        15.0,    // 15°C temperature
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
