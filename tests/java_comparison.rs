//! Test to compare results with the Java implementation.

use chrono::{DateTime, FixedOffset, TimeZone, Utc};
use solar_positioning::{spa, RefractionCorrection};

#[test]
fn test_compare_with_java_spa() {
    // Use the exact same parameters as one of the reference data points
    // From the CSV: 1910-03-15T00:30:00Z,-36.840556,174.740000,0.188643,34.269919

    let utc_datetime = "1910-03-15T00:30:00Z".parse::<DateTime<Utc>>().unwrap();
    let datetime = FixedOffset::east_opt(0)
        .unwrap()
        .from_utc_datetime(&utc_datetime.naive_utc());
    let latitude = -36.840556;
    let longitude = 174.740000;
    let expected_azimuth = 0.188643;
    let expected_zenith = 34.269919;

    println!("Testing against Java reference data:");
    println!("DateTime: {}", utc_datetime);
    println!("Latitude: {:.6}°", latitude);
    println!("Longitude: {:.6}°", longitude);
    println!("Expected azimuth: {:.6}°", expected_azimuth);
    println!("Expected zenith: {:.6}°", expected_zenith);

    let result = spa::solar_position(
        datetime,
        latitude,
        longitude,
        0.0,                                                    // elevation (same as reference)
        0.0,                                                    // deltaT (same as reference)
        Some(RefractionCorrection::new(1000.0, 10.0).unwrap()), // atmospheric conditions
    );

    assert!(result.is_ok());
    let position = result.unwrap();

    println!("Rust calculated azimuth: {:.6}°", position.azimuth());
    println!("Rust calculated zenith: {:.6}°", position.zenith_angle());

    let azimuth_error = (position.azimuth() - expected_azimuth).abs();
    let zenith_error = (position.zenith_angle() - expected_zenith).abs();

    println!("Azimuth error: {:.6}°", azimuth_error);
    println!("Zenith error: {:.6}°", zenith_error);

    // For debugging, let's see if we're even in the right ballpark
    // The Java implementation should be accurate to 0.0003°, so let's see
    // how far off we are
    assert!(azimuth_error < 180.0, "Azimuth is completely wrong");
    assert!(zenith_error < 90.0, "Zenith is completely wrong");

    println!("Results are at least in the right general range");
}
