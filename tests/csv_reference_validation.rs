//! Comprehensive validation against the same CSV reference data used in the Java version.

use chrono::{DateTime, Utc};
use solar_positioning::{grena3, spa};
use std::fs::File;
use std::io::{BufRead, BufReader};

const EPSILON: f64 = 0.001; // Allow 0.001° deviation for floating-point differences

#[test]
fn validate_spa_against_csv_reference_data() {
    let file = File::open("tests/data/test/azimuth_zenith/spa_reference_testdata.csv")
        .expect("SPA reference CSV file should exist");
    let reader = BufReader::new(file);

    let mut test_count = 0;
    let mut max_azimuth_error = 0.0_f64;
    let mut max_zenith_error = 0.0_f64;

    for line in reader.lines() {
        let line = line.unwrap();

        // Skip comments and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        // Parse CSV line: datetime,latitude,longitude,expected_azimuth,expected_zenith
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() != 5 {
            continue;
        }

        let datetime_str = parts[0];
        let latitude: f64 = parts[1].parse().unwrap();
        let longitude: f64 = parts[2].parse().unwrap();
        let expected_azimuth: f64 = parts[3].parse().unwrap();
        let expected_zenith: f64 = parts[4].parse().unwrap();

        let datetime = datetime_str.parse::<DateTime<Utc>>().unwrap();

        // Calculate with our SPA implementation
        // Parameters from CSV header: deltaT=0, pressure=1000mb, temperature=10°C, elevation=0m
        let result = spa::solar_position(
            datetime, latitude, longitude, 0.0,    // elevation
            0.0,    // deltaT
            1000.0, // pressure
            10.0,   // temperature
        );

        assert!(
            result.is_ok(),
            "SPA calculation failed for {}",
            datetime_str
        );
        let position = result.unwrap();

        let azimuth_error = (position.azimuth() - expected_azimuth).abs();
        let zenith_error = (position.zenith_angle() - expected_zenith).abs();

        max_azimuth_error = max_azimuth_error.max(azimuth_error);
        max_zenith_error = max_zenith_error.max(zenith_error);

        // Only print failures to keep output manageable
        if azimuth_error > EPSILON || zenith_error > EPSILON {
            println!(
                "FAIL {}: Az {:.6}° (exp {:.6}°, err {:.6}°), Zen {:.6}° (exp {:.6}°, err {:.6}°)",
                datetime_str,
                position.azimuth(),
                expected_azimuth,
                azimuth_error,
                position.zenith_angle(),
                expected_zenith,
                zenith_error
            );
        }

        assert!(
            azimuth_error < EPSILON,
            "Azimuth error {:.6}° exceeds tolerance {:.6}° for {}",
            azimuth_error,
            EPSILON,
            datetime_str
        );

        assert!(
            zenith_error < EPSILON,
            "Zenith error {:.6}° exceeds tolerance {:.6}° for {}",
            zenith_error,
            EPSILON,
            datetime_str
        );

        test_count += 1;
    }

    println!("✓ Validated {} SPA test cases", test_count);
    println!("✓ Max azimuth error: {:.6}°", max_azimuth_error);
    println!("✓ Max zenith error: {:.6}°", max_zenith_error);
    assert!(test_count > 0, "Should have tested some cases");
}

#[test]
fn validate_grena3_basic_functionality() {
    // Grena3 is designed for 2010-2110, so test with modern dates
    let test_cases = [
        ("2023-06-21T12:00:00Z", 37.7749, -122.4194), // San Francisco, summer solstice
        ("2023-12-21T12:00:00Z", 40.7128, -74.0060),  // New York, winter solstice
        ("2023-03-20T12:00:00Z", 51.5074, -0.1278),   // London, spring equinox
        ("2023-09-23T12:00:00Z", -33.8688, 151.2093), // Sydney, autumn equinox
    ];

    for (datetime_str, latitude, longitude) in test_cases {
        let datetime = datetime_str.parse::<DateTime<Utc>>().unwrap();

        let result = grena3::solar_position(datetime, latitude, longitude, 69.0);

        assert!(
            result.is_ok(),
            "Grena3 calculation failed for {}",
            datetime_str
        );
        let position = result.unwrap();

        // Basic sanity checks
        assert!(position.azimuth() >= 0.0 && position.azimuth() <= 360.0);
        assert!(position.zenith_angle() >= 0.0 && position.zenith_angle() <= 180.0);

        println!(
            "{}: Az {:.3}°, Zen {:.3}°, Elev {:.3}°",
            datetime_str,
            position.azimuth(),
            position.zenith_angle(),
            position.elevation_angle()
        );
    }
}

#[test]
fn compare_spa_vs_grena3_accuracy() {
    // Compare SPA and Grena3 for modern dates where both should work
    let test_cases = [
        ("2023-06-21T12:00:00Z", 37.7749, -122.4194),
        ("2023-12-21T12:00:00Z", 40.7128, -74.0060),
    ];

    for (datetime_str, latitude, longitude) in test_cases {
        let datetime = datetime_str.parse::<DateTime<Utc>>().unwrap();
        let delta_t = 69.0; // 2023 deltaT

        let spa_result = spa::solar_position(
            datetime,
            latitude,
            longitude,
            0.0,
            delta_t,
            f64::NAN,
            f64::NAN,
        )
        .unwrap();
        let grena3_result = grena3::solar_position(datetime, latitude, longitude, delta_t).unwrap();

        let azimuth_diff = (spa_result.azimuth() - grena3_result.azimuth()).abs();
        let zenith_diff = (spa_result.zenith_angle() - grena3_result.zenith_angle()).abs();

        println!(
            "{}: SPA Az {:.4}°/Zen {:.4}°, Grena3 Az {:.4}°/Zen {:.4}° (diff: {:.4}°/{:.4}°)",
            datetime_str,
            spa_result.azimuth(),
            spa_result.zenith_angle(),
            grena3_result.azimuth(),
            grena3_result.zenith_angle(),
            azimuth_diff,
            zenith_diff
        );

        // Grena3 should be within 0.01° (per paper), but allow some margin
        assert!(
            azimuth_diff < 0.02,
            "Azimuth difference too large: {:.6}°",
            azimuth_diff
        );
        assert!(
            zenith_diff < 0.02,
            "Zenith difference too large: {:.6}°",
            zenith_diff
        );
    }
}
