//! Validate SPA implementation against NREL reference data.

use chrono::{DateTime, Utc};
use solar_positioning::spa;

const EPSILON: f64 = 0.001; // Allow 0.001° deviation for floating-point differences

#[test]
fn validate_against_nrel_reference_data() {
    // Reference data from spa_reference_testdata.csv
    // Parameters: deltaT=0, pressure=1000mb, temperature=10°C, elevation=0m
    let test_cases = [
        // Format: (datetime, lat, lon, expected_azimuth, expected_zenith)
        (
            "1910-03-15T00:30:00Z",
            -36.840556,
            174.740000,
            0.188643,
            34.269919,
        ),
        (
            "1910-03-15T03:30:00Z",
            -36.840556,
            174.740000,
            298.894756,
            53.637925,
        ),
        (
            "1910-03-15T06:30:00Z",
            -36.840556,
            174.740000,
            268.082350,
            88.143823,
        ),
        (
            "1910-03-15T09:30:00Z",
            -36.840556,
            174.740000,
            237.156205,
            122.642657,
        ),
        (
            "1910-03-15T12:30:00Z",
            -36.840556,
            174.740000,
            180.112832,
            140.797480,
        ),
    ];

    for (datetime_str, latitude, longitude, expected_azimuth, expected_zenith) in test_cases {
        let datetime = datetime_str.parse::<DateTime<Utc>>().unwrap();

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

        println!(
            "{}: Azimuth {:.6}° (expected {:.6}°, error {:.6}°), Zenith {:.6}° (expected {:.6}°, error {:.6}°)",
            datetime_str,
            position.azimuth(),
            expected_azimuth,
            azimuth_error,
            position.zenith_angle(),
            expected_zenith,
            zenith_error
        );

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
    }
}

#[test]
fn test_modern_date_calculation() {
    // Test a modern date to ensure our fix works
    // Use solar noon time for San Francisco (approximately 20:00 UTC during summer solstice)
    let datetime = "2023-06-21T20:00:00Z".parse::<DateTime<Utc>>().unwrap();

    let result = spa::solar_position(
        datetime, 37.7749,   // San Francisco latitude
        -122.4194, // San Francisco longitude
        0.0,       // elevation
        69.0,      // deltaT for 2023
        1013.25,   // standard pressure
        15.0,      // temperature
    );

    assert!(
        result.is_ok(),
        "SPA calculation should work for modern dates"
    );
    let position = result.unwrap();

    // At summer solstice noon in San Francisco, sun should be:
    // - Nearly due south (azimuth ~180°)
    // - At high elevation (low zenith angle, around 20-30°)

    println!(
        "2023-06-21T20:00:00Z San Francisco: Azimuth {:.6}°, Zenith {:.6}°, Elevation {:.6}°",
        position.azimuth(),
        position.zenith_angle(),
        position.elevation_angle()
    );

    // Sanity checks for summer solstice at noon
    assert!(
        position.azimuth() > 150.0 && position.azimuth() < 210.0,
        "Azimuth should be roughly south at solar noon"
    );
    assert!(
        position.zenith_angle() > 0.0 && position.zenith_angle() < 45.0,
        "Zenith angle should be small at summer solstice noon"
    );
    assert!(
        position.elevation_angle() > 45.0 && position.elevation_angle() < 90.0,
        "Elevation should be high at summer solstice noon"
    );
}
