//! Test the standard sunce performance case: 50:55:0.1 10:15:0.1 2024 position --step=3h

#![cfg(feature = "unstable")]

use chrono::{DateTime, Duration, Utc};
use solar_positioning::spa;
use std::collections::HashMap;
use std::time::Instant;

#[test]
fn test_sunce_standard_case_optimization() {
    // sunce --perf --format=CSV --no-headers 50:55:0.1 10:15:0.1 2024 position --step=3h

    // Generate coordinate grid: 50:55:0.1 (lat) × 10:15:0.1 (lon)
    let mut coordinates = Vec::new();
    let mut lat = 50.0;
    while lat <= 55.0 + 1e-9 {
        // Add epsilon for floating point precision
        let mut lon = 10.0;
        while lon <= 15.0 + 1e-9 {
            coordinates.push((lat, lon));
            lon += 0.1;
        }
        lat += 0.1;
    }

    // Generate time series: 2024 with --step=3h
    let start_time = "2024-01-01T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let end_time = "2025-01-01T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let mut times = Vec::new();
    let mut current_time = start_time;
    while current_time < end_time {
        times.push(current_time);
        current_time += Duration::hours(3);
    }

    let coord_count = coordinates.len();
    let time_count = times.len();
    let total_calculations = coord_count * time_count;

    println!("Sunce standard case:");
    println!(
        "Coordinates: {} ({}×{})",
        coord_count,
        ((55.0 - 50.0) / 0.1 + 1.0) as i32,
        ((15.0 - 10.0) / 0.1 + 1.0) as i32
    );
    println!("Time points: {} (3-hour steps through 2024)", time_count);
    println!(
        "Total calculations: {:.1}M",
        total_calculations as f64 / 1_000_000.0
    );

    // Naive approach: full SPA for each coordinate+time pair
    println!("\n--- Naive Approach ---");
    let start = Instant::now();
    let mut naive_count = 0;

    // Only test a subset to keep test runtime reasonable
    let test_coords = &coordinates[0..100]; // 100 coordinates
    let test_times = &times[0..100]; // 100 time points = 10,000 calculations

    for &time in test_times {
        for &(lat, lon) in test_coords {
            let _result = spa::solar_position(
                time,
                lat,
                lon,
                0.0,
                69.0,
                Some(solar_positioning::RefractionCorrection::new(1013.25, 15.0).unwrap()),
            )
            .unwrap();
            naive_count += 1;
        }
    }
    let naive_duration = start.elapsed();

    // Time-cached approach
    println!("--- Time-Cached Approach ---");
    let start = Instant::now();
    let mut cached_count = 0;
    let mut time_cache: HashMap<DateTime<Utc>, spa::SpaTimeDependent> = HashMap::new();

    for &time in test_times {
        let time_parts = time_cache
            .entry(time)
            .or_insert_with(|| spa::spa_time_dependent_parts(time, 69.0).unwrap());

        for &(lat, lon) in test_coords {
            let _result = spa::spa_with_time_dependent_parts(
                lat,
                lon,
                0.0,
                Some(solar_positioning::RefractionCorrection::new(1013.25, 15.0).unwrap()),
                time_parts,
            )
            .unwrap();
            cached_count += 1;
        }
    }
    let cached_duration = start.elapsed();

    assert_eq!(naive_count, cached_count);
    assert_eq!(naive_count, test_coords.len() * test_times.len());

    let speedup = naive_duration.as_secs_f64() / cached_duration.as_secs_f64();

    println!("Test subset: {} calculations", naive_count);
    println!("Naive: {:.3}s", naive_duration.as_secs_f64());
    println!("Cached: {:.3}s", cached_duration.as_secs_f64());
    println!("Speedup: {:.2}x", speedup);
    println!("Cache entries: {}", time_cache.len());

    // Project to full sunce case
    let full_naive_time =
        naive_duration.as_secs_f64() * (total_calculations as f64 / naive_count as f64);
    let full_cached_time =
        cached_duration.as_secs_f64() * (total_calculations as f64 / naive_count as f64);

    println!(
        "\n--- Projected Full Case ({:.1}M calculations) ---",
        total_calculations as f64 / 1_000_000.0
    );
    println!("Naive projection: {:.1}s", full_naive_time);
    println!("Cached projection: {:.1}s", full_cached_time);
    println!(
        "Projected speedup: {:.2}x",
        full_naive_time / full_cached_time
    );

    // Verify we get meaningful speedup
    assert!(
        speedup > 1.5,
        "Expected speedup > 1.5x, got {:.2}x",
        speedup
    );
}

#[test]
fn test_coordinate_grid_generation() {
    // Verify we generate the coordinate grid correctly (matches sunce behavior)
    let mut coordinates = Vec::new();
    let mut lat = 50.0;
    while lat <= 55.0 + 1e-9 {
        let mut lon = 10.0;
        while lon <= 15.0 + 1e-9 {
            coordinates.push((lat, lon));
            lon += 0.1;
        }
        lat += 0.1;
    }

    // Should be 51×51 = 2,601 coordinates
    let expected_lat_count = ((55.0 - 50.0) / 0.1 + 1.0) as usize;
    let expected_lon_count = ((15.0 - 10.0) / 0.1 + 1.0) as usize;
    let expected_total = expected_lat_count * expected_lon_count;

    assert_eq!(coordinates.len(), expected_total);
    assert_eq!(coordinates.len(), 2601); // 51×51

    // Check corners (allow small floating point errors)
    assert_eq!(coordinates[0], (50.0, 10.0));
    let last = coordinates[coordinates.len() - 1];
    assert!((last.0 - 55.0_f64).abs() < 1e-10);
    assert!((last.1 - 15.0_f64).abs() < 1e-10);
}

#[test]
fn test_time_series_generation() {
    // Verify we generate the time series correctly (matches sunce 2024 --step=3h)
    let start_time = "2024-01-01T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let end_time = "2025-01-01T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let mut times = Vec::new();
    let mut current_time = start_time;
    while current_time < end_time {
        times.push(current_time);
        current_time += Duration::hours(3);
    }

    // 2024 is a leap year: 366 days × 8 times per day = 2,928 time points
    let expected_count = 366 * 8;
    assert_eq!(times.len(), expected_count);

    // Check first and last
    assert_eq!(times[0], start_time);
    assert_eq!(
        times[times.len() - 1],
        "2024-12-31T21:00:00Z".parse::<DateTime<Utc>>().unwrap()
    );
}
