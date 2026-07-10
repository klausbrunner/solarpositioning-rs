//! Test to compare results with the Java implementation.

#![cfg(feature = "chrono")]

use chrono::{DateTime, FixedOffset, TimeZone, Utc};
use solar_positioning::{spa, RefractionCorrection};

const TOLERANCE_DEGREES: f64 = 0.0001;

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

    let position = spa::solar_position(
        datetime,
        latitude,
        longitude,
        0.0,                                                    // elevation (same as reference)
        0.0,                                                    // deltaT (same as reference)
        Some(RefractionCorrection::new(1000.0, 10.0).unwrap()), // atmospheric conditions
    )
    .unwrap();

    let azimuth_difference = (position.azimuth() - expected_azimuth).abs();
    let azimuth_error = azimuth_difference.min(360.0 - azimuth_difference);
    let zenith_error = (position.zenith_angle() - expected_zenith).abs();

    assert!(
        azimuth_error < TOLERANCE_DEGREES,
        "azimuth error {azimuth_error:.6}° exceeds {TOLERANCE_DEGREES:.6}°"
    );
    assert!(
        zenith_error < TOLERANCE_DEGREES,
        "zenith error {zenith_error:.6}° exceeds {TOLERANCE_DEGREES:.6}°"
    );
}
