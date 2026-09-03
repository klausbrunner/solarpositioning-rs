//! Validate SPA positions against NREL reference data.

#![cfg(all(feature = "chrono", feature = "std"))]

use chrono::{DateTime, Utc};
use csv::ReaderBuilder;
use solar_positioning::{spa, RefractionCorrection};
use std::error::Error;
use std::fs::File;

const TOLERANCE_DEGREES: f64 = 0.0003;

#[test]
fn spa_matches_nrel_reference_data() -> Result<(), Box<dyn Error>> {
    let file = File::open("tests/data/spa_reference_testdata.csv")?;
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(file);
    let refraction = RefractionCorrection::new(1000.0, 10.0)?;

    let mut count = 0;
    for record in reader.records() {
        let record = record?;
        let datetime = record[0].parse::<DateTime<Utc>>()?;
        let latitude = record[1].parse()?;
        let longitude = record[2].parse()?;
        let expected_azimuth: f64 = record[3].parse()?;
        let expected_zenith: f64 = record[4].parse()?;

        let position =
            spa::solar_position(datetime, latitude, longitude, 0.0, 0.0, Some(refraction))?;
        let azimuth_error = (position.azimuth() - expected_azimuth).abs();
        let zenith_error = (position.zenith_angle() - expected_zenith).abs();

        assert!(
            azimuth_error < TOLERANCE_DEGREES,
            "azimuth error {azimuth_error:.6}° for {datetime}"
        );
        assert!(
            zenith_error < TOLERANCE_DEGREES,
            "zenith error {zenith_error:.6}° for {datetime}"
        );
        count += 1;
    }

    assert!(count > 0, "no reference records loaded");
    Ok(())
}
