//! Test sunrise/sunset calculations against reference data.

use chrono::{DateTime, Datelike, FixedOffset, TimeZone, Timelike, Utc};
use csv::ReaderBuilder;
use solar_positioning::{spa, types::SunriseResult, RefractionCorrection};
use std::error::Error;
use std::fs::File;

#[derive(Debug)]
struct SunriseTestRecord {
    datetime: DateTime<Utc>,
    latitude: f64,
    longitude: f64,
    expected_sunset: String,
    expected_transit: String,
    expected_sunrise: String,
}

impl SunriseTestRecord {
    fn from_csv_record(record: &csv::StringRecord) -> Result<Self, Box<dyn Error>> {
        Ok(Self {
            datetime: record[0].parse()?,
            latitude: record[1].parse()?,
            longitude: record[2].parse()?,
            expected_sunrise: record[3].to_string(),
            expected_transit: record[4].to_string(),
            expected_sunset: record[5].to_string(),
        })
    }
}

fn time_difference_seconds(expected: &str, actual: &DateTime<Utc>) -> f64 {
    // Parse expected time in HH:MM:SS format
    let parts: Vec<&str> = expected.split(':').collect();
    if parts.len() != 3 {
        return f64::INFINITY;
    }

    let expected_hours: u32 = parts[0].parse().unwrap_or(0);
    let expected_minutes: u32 = parts[1].parse().unwrap_or(0);
    let expected_seconds: u32 = parts[2].parse().unwrap_or(0);

    let expected_total_seconds = expected_hours * 3600 + expected_minutes * 60 + expected_seconds;
    let actual_total_seconds = actual.hour() * 3600 + actual.minute() * 60 + actual.second();

    (actual_total_seconds as i64 - expected_total_seconds as i64).abs() as f64
}

#[test]
fn test_sunrise_sunset_debug_single_case() -> Result<(), Box<dyn Error>> {
    // Debug a single test case
    let test_datetime_utc = "1910-03-05T12:30:00Z".parse::<DateTime<Utc>>()?;
    let test_datetime = FixedOffset::east_opt(0)
        .unwrap()
        .from_utc_datetime(&test_datetime_utc.naive_utc());
    let test_latitude = -36.840556;
    let test_longitude = 174.740000;
    let delta_t = 0.0;

    println!("=== DEBUG SINGLE CASE ===");
    println!("Input: {}", test_datetime_utc);
    println!("Lat: {:.6}, Lon: {:.6}", test_latitude, test_longitude);

    // Java SPA returns: Sunrise=18:09:57, Transit=00:32:56, Sunset=06:56:15
    // CSV format: date lat lon sunrise transit sunset
    // CSV data: 1910-03-05T12:30:00Z,-36.840556,174.740000,18:09:56,00:32:56,06:56:15
    let expected_sunrise = "18:09:56";
    let expected_transit = "00:32:56";
    let expected_sunset = "06:56:15";

    let result = spa::sunrise_sunset(
        test_datetime,
        test_latitude,
        test_longitude,
        delta_t,
        -0.833, // Standard sunrise/sunset elevation angle
    )?;

    match result {
        SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        } => {
            println!(
                "Expected: Sunrise={}, Transit={}, Sunset={}",
                expected_sunrise, expected_transit, expected_sunset
            );
            println!(
                "Actual:   Sunrise={}, Transit={}, Sunset={}",
                sunrise.format("%H:%M:%S"),
                transit.format("%H:%M:%S"),
                sunset.format("%H:%M:%S")
            );

            // Show the Julian date calculation details
            let day_start = test_datetime_utc
                .date_naive()
                .and_hms_opt(0, 0, 0)
                .unwrap()
                .and_utc();
            println!("Day start: {}", day_start);

            // Let's also check what day_start gives us
            let day_start_fixed = FixedOffset::east_opt(0)
                .unwrap()
                .from_utc_datetime(&day_start.naive_utc());
            let jd = solar_positioning::time::JulianDate::from_datetime(&day_start_fixed, delta_t)?;
            println!("Julian Date: {:.6}", jd.julian_date());
            println!("JME: {:.6}", jd.julian_ephemeris_millennium());

            // Let's check if our basic SPA position calculation works for the same inputs
            let position = spa::solar_position(
                test_datetime,
                test_latitude,
                test_longitude,
                0.0, // elevation = 0
                delta_t,
                Some(RefractionCorrection::standard()), // atmospheric conditions
            )?;
            println!(
                "SPA Position at input time: Azimuth={:.3}°, Zenith={:.3}°",
                position.azimuth(),
                position.zenith_angle()
            );

            // Let's check position at the calculated sunrise time
            let sunrise_position = spa::solar_position(
                sunrise,
                test_latitude,
                test_longitude,
                0.0,
                delta_t,
                Some(RefractionCorrection::standard()),
            )?;
            println!(
                "SPA Position at calculated sunrise: Azimuth={:.3}°, Elevation={:.3}°",
                sunrise_position.azimuth(),
                sunrise_position.elevation_angle()
            );

            // At sunrise, elevation should be around -0.833 degrees
            println!("Expected elevation at sunrise: -0.833°");

            // Let's also check our position calculation at the expected transit time
            let expected_transit_time_utc = "1910-03-05T00:32:56Z".parse::<DateTime<Utc>>()?;
            let expected_transit_time = FixedOffset::east_opt(0)
                .unwrap()
                .from_utc_datetime(&expected_transit_time_utc.naive_utc());
            let transit_position = spa::solar_position(
                expected_transit_time,
                test_latitude,
                test_longitude,
                0.0,
                delta_t,
                Some(RefractionCorrection::standard()),
            )?;
            println!(
                "Expected transit time {} -> Azimuth={:.5}°, Zenith={:.5}°",
                expected_transit_time.format("%H:%M:%S"),
                transit_position.azimuth(),
                transit_position.zenith_angle()
            );
            println!("Java reference at transit: Azimuth=0.00364°, Zenith=30.36656°");
        }
        _ => {
            panic!("Expected regular day result");
        }
    }

    Ok(())
}

#[test]
fn test_sunrise_sunset_against_spa_reference_data() -> Result<(), Box<dyn Error>> {
    let file = File::open("tests/data/test/sunrise/spa_reference_testdata.csv")?;
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(file);

    let mut records = Vec::new();
    for result in reader.records() {
        let record = result?;
        if !record.is_empty() && record.len() >= 6 {
            match SunriseTestRecord::from_csv_record(&record) {
                Ok(test_record) => records.push(test_record),
                Err(e) => {
                    eprintln!("Skipping invalid record: {:?}, error: {}", record, e);
                    continue;
                }
            }
        }
    }

    println!("Loaded {} sunrise/sunset test records", records.len());

    let mut max_sunrise_error = 0.0f64;
    let mut max_sunset_error = 0.0f64;
    let mut max_transit_error = 0.0f64;
    let mut failed_cases = 0;

    for (i, record) in records.iter().enumerate() {
        // Use Delta T = 0 to match reference data
        let delta_t = 0.0;

        let record_datetime_fixed = FixedOffset::east_opt(0)
            .unwrap()
            .from_utc_datetime(&record.datetime.naive_utc());
        match spa::sunrise_sunset(
            record_datetime_fixed,
            record.latitude,
            record.longitude,
            delta_t,
            -0.833, // Standard sunrise/sunset elevation angle
        ) {
            Ok(result) => {
                match result {
                    SunriseResult::RegularDay {
                        sunrise,
                        transit,
                        sunset,
                    } => {
                        let sunrise_utc = sunrise.with_timezone(&Utc);
                        let sunrise_error =
                            time_difference_seconds(&record.expected_sunrise, &sunrise_utc);
                        let transit_utc = transit.with_timezone(&Utc);
                        let transit_error =
                            time_difference_seconds(&record.expected_transit, &transit_utc);
                        let sunset_utc = sunset.with_timezone(&Utc);
                        let sunset_error =
                            time_difference_seconds(&record.expected_sunset, &sunset_utc);

                        max_sunrise_error = max_sunrise_error.max(sunrise_error);
                        max_transit_error = max_transit_error.max(transit_error);
                        max_sunset_error = max_sunset_error.max(sunset_error);

                        // Different tolerances for different calculations:
                        // Transit times should be most accurate (sun at meridian, minimal atmospheric effects)
                        let transit_tolerance = 1.0; // 1 second for transit
                                                     // Sunrise/sunset have more atmospheric uncertainty at horizon
                        let horizon_tolerance = 120.0; // 2 minutes for sunrise/sunset

                        if sunrise_error > horizon_tolerance
                            || sunset_error > horizon_tolerance
                            || transit_error > transit_tolerance
                        {
                            println!(
                                "Record {}: DateTime={}, Lat={:.6}, Lon={:.6}",
                                i + 1,
                                record.datetime,
                                record.latitude,
                                record.longitude
                            );
                            println!(
                                "  Expected: Sunrise={}, Transit={}, Sunset={}",
                                record.expected_sunrise,
                                record.expected_transit,
                                record.expected_sunset
                            );
                            println!(
                                "  Actual:   Sunrise={}, Transit={}, Sunset={}",
                                sunrise.format("%H:%M:%S"),
                                transit.format("%H:%M:%S"),
                                sunset.format("%H:%M:%S")
                            );
                            println!(
                                "  Errors:   Sunrise={:.1}s, Transit={:.1}s, Sunset={:.1}s",
                                sunrise_error, transit_error, sunset_error
                            );
                            failed_cases += 1;
                        }
                    }
                    _ => {
                        // Skip polar day/night cases for now as reference data only contains regular days
                        continue;
                    }
                }
            }
            Err(e) => {
                println!(
                    "SPA sunrise/sunset calculation failed for record {}: {}",
                    i + 1,
                    e
                );
                failed_cases += 1;
            }
        }
    }

    println!("Maximum sunrise error: {:.1} seconds", max_sunrise_error);
    println!("Maximum transit error: {:.1} seconds", max_transit_error);
    println!("Maximum sunset error: {:.1} seconds", max_sunset_error);
    println!("Failed cases: {}", failed_cases);

    // Different tolerances for different calculations
    let transit_tolerance = 1.0; // Transit should be most accurate
    let horizon_tolerance = 120.0; // Sunrise/sunset have atmospheric uncertainty

    assert!(
        max_sunrise_error < horizon_tolerance,
        "Maximum sunrise error {:.1}s exceeds tolerance {:.1}s",
        max_sunrise_error,
        horizon_tolerance
    );

    assert!(
        max_transit_error < transit_tolerance,
        "Maximum transit error {:.1}s exceeds tolerance {:.1}s",
        max_transit_error,
        transit_tolerance
    );

    assert!(
        max_sunset_error < horizon_tolerance,
        "Maximum sunset error {:.1}s exceeds tolerance {:.1}s",
        max_sunset_error,
        horizon_tolerance
    );

    Ok(())
}

#[test]
fn test_polar_transit_accuracy_svalbard() -> Result<(), Box<dyn Error>> {
    // Regression test for polar latitude transit accuracy
    // Svalbard, Norway (78.0°N, 15.0°E) on 2024-01-01
    // Verify chrono wrapper stays aligned with SPA's 0 UT anchored algorithm.
    let latitude = 78.0;
    let longitude = 15.0;
    let tz = FixedOffset::east_opt(3600).unwrap();
    let date = tz.with_ymd_and_hms(2024, 1, 1, 12, 0, 0).unwrap();

    let result = spa::sunrise_sunset(date, latitude, longitude, 0.0, -0.833)?;

    match result {
        SunriseResult::AllNight { transit } => {
            // Chrono wrapper treats the local calendar date as the SPA UTC calculation date.
            let base_utc_date = chrono::NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
            let expected = match spa::sunrise_sunset_utc(
                base_utc_date.year(),
                base_utc_date.month(),
                base_utc_date.day(),
                latitude,
                longitude,
                0.0,
                -0.833,
            )? {
                SunriseResult::AllNight { transit } => transit,
                _ => panic!("expected polar night from UTC API"),
            };
            let expected_utc = base_utc_date.and_hms_opt(0, 0, 0).unwrap().and_utc()
                + chrono::Duration::milliseconds((expected.hours() * 3_600_000.0) as i64);

            assert_eq!(
                transit.with_timezone(&Utc).naive_utc(),
                expected_utc.naive_utc()
            );

            // Verify azimuth is very close to 180° at transit
            let pos = spa::solar_position(
                transit,
                latitude,
                longitude,
                0.0,
                0.0,
                Some(RefractionCorrection::standard()),
            )?;
            let azimuth_error = (pos.azimuth() - 180.0).abs();
            assert!(
                azimuth_error < 0.001,
                "Azimuth error {:.6}° exceeds tolerance",
                azimuth_error
            );
        }
        _ => panic!(
            "Expected polar night for Svalbard in January, got: {:?}",
            result
        ),
    }

    Ok(())
}
