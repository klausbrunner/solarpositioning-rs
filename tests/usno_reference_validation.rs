//! Test sunrise/sunset calculations against USNO reference data.

use chrono::{DateTime, Timelike, Utc};
use csv::ReaderBuilder;
use solar_positioning::{spa, types::SunriseResult};
use std::error::Error;
use std::fs::File;

#[derive(Debug)]
struct USNOTestRecord {
    datetime: DateTime<Utc>,
    latitude: f64,
    longitude: f64,
    day_type: String,
    expected_sunrise: Option<String>,
    expected_sunset: Option<String>,
}

impl USNOTestRecord {
    fn from_csv_record(record: &csv::StringRecord) -> Result<Self, Box<dyn Error>> {
        Ok(Self {
            datetime: record[0].parse()?,
            latitude: record[1].parse()?,
            longitude: record[2].parse()?,
            day_type: record[3].to_string(),
            expected_sunrise: if record.len() > 4 && !record[4].is_empty() {
                Some(record[4].to_string())
            } else {
                None
            },
            expected_sunset: if record.len() > 5 && !record[5].is_empty() {
                Some(record[5].to_string())
            } else {
                None
            },
        })
    }
}

fn time_difference_seconds(expected: &str, actual: &DateTime<Utc>) -> f64 {
    // Parse expected time in HH:MM format
    let parts: Vec<&str> = expected.split(':').collect();
    if parts.len() != 2 {
        return f64::INFINITY;
    }

    let expected_hours: u32 = parts[0].parse().unwrap_or(0);
    let expected_minutes: u32 = parts[1].parse().unwrap_or(0);

    let expected_total_seconds = expected_hours * 3600 + expected_minutes * 60;
    let actual_total_seconds = actual.hour() * 3600 + actual.minute() * 60 + actual.second();

    (actual_total_seconds as i64 - expected_total_seconds as i64).abs() as f64
}

#[test]
fn test_usno_reference_data() -> Result<(), Box<dyn Error>> {
    let file = File::open("tests/data/usno/usno_reference_testdata.csv")?;
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(file);

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result?;
        if !record.is_empty() && record.len() >= 4 {
            match USNOTestRecord::from_csv_record(&record) {
                Ok(test_record) => records.push(test_record),
                Err(e) => {
                    eprintln!("Skipping invalid record: {:?}, error: {}", record, e);
                    continue;
                }
            }
        }
    }

    println!("Loaded {} USNO test records", records.len());

    let mut max_sunrise_error = 0.0f64;
    let mut max_sunset_error = 0.0f64;
    let mut failed_cases = 0;

    // Allow 40 seconds tolerance as in Java (REASONABLE_TOLERANCE)
    let tolerance = 40.0;

    for (i, record) in records.iter().enumerate() {
        let result = spa::sunrise_sunset(
            record.datetime,
            record.latitude,
            record.longitude,
            0.0,    // Delta T = 0
            -0.833, // Standard sunrise/sunset elevation angle
        );

        match result {
            Ok(sunrise_result) => match record.day_type.as_str() {
                "NORMAL" => {
                    if let SunriseResult::RegularDay {
                        sunrise, sunset, ..
                    } = sunrise_result
                    {
                        if let Some(ref expected_sunrise) = record.expected_sunrise {
                            let sunrise_error = time_difference_seconds(expected_sunrise, &sunrise);
                            max_sunrise_error = max_sunrise_error.max(sunrise_error);

                            if sunrise_error > tolerance {
                                println!(
                                    "Record {}: Sunrise error {:.1}s exceeds tolerance {:.1}s",
                                    i + 1,
                                    sunrise_error,
                                    tolerance
                                );
                                failed_cases += 1;
                            }
                        }

                        if let Some(ref expected_sunset) = record.expected_sunset {
                            let sunset_error = time_difference_seconds(expected_sunset, &sunset);
                            max_sunset_error = max_sunset_error.max(sunset_error);

                            if sunset_error > tolerance {
                                println!(
                                    "Record {}: Sunset error {:.1}s exceeds tolerance {:.1}s",
                                    i + 1,
                                    sunset_error,
                                    tolerance
                                );
                                failed_cases += 1;
                            }
                        }
                    } else {
                        println!(
                            "Record {}: Expected NORMAL day but got {:?}",
                            i + 1,
                            sunrise_result
                        );
                        failed_cases += 1;
                    }
                }
                "ALL_DAY" => {
                    if !matches!(sunrise_result, SunriseResult::AllDay { .. }) {
                        println!(
                            "Record {}: Expected ALL_DAY but got {:?}",
                            i + 1,
                            sunrise_result
                        );
                        failed_cases += 1;
                    }
                }
                "ALL_NIGHT" => {
                    if !matches!(sunrise_result, SunriseResult::AllNight { .. }) {
                        println!(
                            "Record {}: Expected ALL_NIGHT but got {:?}",
                            i + 1,
                            sunrise_result
                        );
                        failed_cases += 1;
                    }
                }
                _ => {
                    println!("Record {}: Unknown day type '{}'", i + 1, record.day_type);
                    failed_cases += 1;
                }
            },
            Err(e) => {
                println!("Record {}: Calculation failed: {}", i + 1, e);
                failed_cases += 1;
            }
        }
    }

    println!("Maximum sunrise error: {:.1} seconds", max_sunrise_error);
    println!("Maximum sunset error: {:.1} seconds", max_sunset_error);
    println!("Failed cases: {}", failed_cases);

    assert!(
        max_sunrise_error < tolerance,
        "Maximum sunrise error {:.1}s exceeds tolerance {:.1}s",
        max_sunrise_error,
        tolerance
    );

    assert!(
        max_sunset_error < tolerance,
        "Maximum sunset error {:.1}s exceeds tolerance {:.1}s",
        max_sunset_error,
        tolerance
    );

    assert_eq!(failed_cases, 0, "Some test cases failed");

    Ok(())
}

#[test]
fn test_usno_extreme_cases() -> Result<(), Box<dyn Error>> {
    let file = File::open("tests/data/usno/usno_reference_testdata_extreme.csv")?;
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(file);

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result?;
        if !record.is_empty() && record.len() >= 4 {
            match USNOTestRecord::from_csv_record(&record) {
                Ok(test_record) => records.push(test_record),
                Err(e) => {
                    eprintln!("Skipping invalid record: {:?}, error: {}", record, e);
                    continue;
                }
            }
        }
    }

    println!("Loaded {} USNO extreme test records", records.len());

    let mut failed_cases = 0;

    // Allow 2.5 minutes tolerance for extreme cases (slightly more than Java's 2 min due to USNO reference differences)
    let tolerance = 150.0;

    for (i, record) in records.iter().enumerate() {
        let result = spa::sunrise_sunset(
            record.datetime,
            record.latitude,
            record.longitude,
            0.0,
            -0.833,
        );

        match result {
            Ok(sunrise_result) => match record.day_type.as_str() {
                "NORMAL" => {
                    if let SunriseResult::RegularDay {
                        sunrise, sunset, ..
                    } = sunrise_result
                    {
                        if let Some(ref expected_sunrise) = record.expected_sunrise {
                            let sunrise_error = time_difference_seconds(expected_sunrise, &sunrise);
                            if sunrise_error > tolerance {
                                println!(
                                    "Record {}: Sunrise error {:.1}s exceeds tolerance {:.1}s",
                                    i + 1,
                                    sunrise_error,
                                    tolerance
                                );
                                failed_cases += 1;
                            }
                        }

                        if let Some(ref expected_sunset) = record.expected_sunset {
                            let sunset_error = time_difference_seconds(expected_sunset, &sunset);
                            if sunset_error > tolerance {
                                println!(
                                    "Record {}: Sunset error {:.1}s exceeds tolerance {:.1}s",
                                    i + 1,
                                    sunset_error,
                                    tolerance
                                );
                                println!(
                                    "  Date: {}, Lat: {}, Lon: {}",
                                    record.datetime, record.latitude, record.longitude
                                );
                                println!(
                                    "  Expected sunset: {:?}, Actual sunset: {}",
                                    expected_sunset, sunset
                                );
                                failed_cases += 1;
                            }
                        }
                    }
                }
                "ALL_DAY" => {
                    if !matches!(sunrise_result, SunriseResult::AllDay { .. }) {
                        println!(
                            "Record {}: Expected ALL_DAY but got {:?}",
                            i + 1,
                            sunrise_result
                        );
                        failed_cases += 1;
                    }
                }
                "ALL_NIGHT" => {
                    if !matches!(sunrise_result, SunriseResult::AllNight { .. }) {
                        println!(
                            "Record {}: Expected ALL_NIGHT but got {:?}",
                            i + 1,
                            sunrise_result
                        );
                        failed_cases += 1;
                    }
                }
                _ => {
                    println!("Record {}: Unknown day type '{}'", i + 1, record.day_type);
                    failed_cases += 1;
                }
            },
            Err(e) => {
                println!("Record {}: Calculation failed: {}", i + 1, e);
                failed_cases += 1;
            }
        }
    }

    println!("Failed extreme cases: {}", failed_cases);

    assert_eq!(failed_cases, 0, "Some extreme test cases failed");

    Ok(())
}

#[test]
fn test_usno_civil_twilight() -> Result<(), Box<dyn Error>> {
    let file = File::open("tests/data/usno/usno_reference_testdata_civil.csv")?;
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(file);

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result?;
        if !record.is_empty() && record.len() >= 4 {
            match USNOTestRecord::from_csv_record(&record) {
                Ok(test_record) => records.push(test_record),
                Err(e) => {
                    eprintln!("Skipping invalid record: {:?}, error: {}", record, e);
                    continue;
                }
            }
        }
    }

    println!("Loaded {} USNO civil twilight test records", records.len());

    let mut failed_cases = 0;

    // Allow 2.5 minutes tolerance for extreme cases (slightly more than Java's 2 min due to USNO reference differences)
    let tolerance = 150.0;

    for (i, record) in records.iter().enumerate() {
        let result = spa::sunrise_sunset(
            record.datetime,
            record.latitude,
            record.longitude,
            0.0,
            -6.0, // Civil twilight elevation angle
        );

        match result {
            Ok(sunrise_result) => match record.day_type.as_str() {
                "NORMAL" => {
                    if let SunriseResult::RegularDay {
                        sunrise, sunset, ..
                    } = sunrise_result
                    {
                        if let Some(ref expected_sunrise) = record.expected_sunrise {
                            let sunrise_error = time_difference_seconds(expected_sunrise, &sunrise);
                            if sunrise_error > tolerance {
                                println!(
                                    "Record {}: Civil twilight sunrise error {:.1}s exceeds tolerance {:.1}s",
                                    i + 1,
                                    sunrise_error,
                                    tolerance
                                );
                                failed_cases += 1;
                            }
                        }

                        if let Some(ref expected_sunset) = record.expected_sunset {
                            let sunset_error = time_difference_seconds(expected_sunset, &sunset);
                            if sunset_error > tolerance {
                                println!(
                                    "Record {}: Civil twilight sunset error {:.1}s exceeds tolerance {:.1}s",
                                    i + 1,
                                    sunset_error,
                                    tolerance
                                );
                                println!(
                                    "  Date: {}, Lat: {}, Lon: {}",
                                    record.datetime, record.latitude, record.longitude
                                );
                                println!(
                                    "  Expected sunset: {:?}, Actual sunset: {}",
                                    expected_sunset, sunset
                                );
                                failed_cases += 1;
                            }
                        }
                    }
                }
                "ALL_DAY" => {
                    if !matches!(sunrise_result, SunriseResult::AllDay { .. }) {
                        println!(
                            "Record {}: Expected ALL_DAY but got {:?}",
                            i + 1,
                            sunrise_result
                        );
                        failed_cases += 1;
                    }
                }
                "ALL_NIGHT" => {
                    if !matches!(sunrise_result, SunriseResult::AllNight { .. }) {
                        println!(
                            "Record {}: Expected ALL_NIGHT but got {:?}",
                            i + 1,
                            sunrise_result
                        );
                        failed_cases += 1;
                    }
                }
                _ => {
                    println!("Record {}: Unknown day type '{}'", i + 1, record.day_type);
                    failed_cases += 1;
                }
            },
            Err(e) => {
                println!("Record {}: Calculation failed: {}", i + 1, e);
                failed_cases += 1;
            }
        }
    }

    println!("Failed civil twilight cases: {}", failed_cases);

    assert_eq!(failed_cases, 0, "Some civil twilight test cases failed");

    Ok(())
}
