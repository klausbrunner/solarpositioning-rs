//! Full test of sunce standard case - actual measurements, both approaches

#![cfg(feature = "unstable")]

use chrono::{DateTime, Duration, Utc};
use solar_positioning::spa;
use std::collections::HashMap;
use std::time::Instant;

#[test]
#[ignore] // Use --ignored to run this expensive test
fn test_sunce_full_case() {
    // Generate full coordinate grid: 50:55:0.1 (lat) × 10:15:0.1 (lon)
    let mut coordinates = Vec::new();
    let mut lat = 50.0;
    while lat <= 55.0 + 1e-9 {
        // Add tiny epsilon for floating point precision
        let mut lon = 10.0;
        while lon <= 15.0 + 1e-9 {
            coordinates.push((lat, lon));
            lon += 0.1;
        }
        lat += 0.1;
    }

    // Generate FULL 2024 with --step=3h (2928 time points = 366 days × 8 per day)
    let start_time = "2024-01-01T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let end_time = "2025-01-01T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let mut times = Vec::new();
    let mut current_time = start_time;
    while current_time < end_time {
        times.push(current_time);
        current_time += Duration::hours(3);
    }

    let total_calculations = coordinates.len() * times.len();
    println!("=== Full Sunce Test Case ===");
    println!(
        "Coordinates: {} (51×51 grid from 50:55:0.1, 10:15:0.1)",
        coordinates.len()
    );
    println!("Time points: {} (full 2024, 3-hour steps)", times.len());
    println!(
        "Total calculations: {} ({:.1}M)",
        total_calculations,
        total_calculations as f64 / 1_000_000.0
    );

    // Naive approach: full SPA for each coordinate+time pair
    println!("\n--- Naive Approach ---");
    let start = Instant::now();
    let mut naive_count = 0;

    for &time in &times {
        for &(lat, lon) in &coordinates {
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

    println!(
        "Completed {} calculations in {:.3}s",
        naive_count,
        naive_duration.as_secs_f64()
    );
    println!(
        "Rate: {:.0} calculations/second",
        naive_count as f64 / naive_duration.as_secs_f64()
    );

    // Time-cached approach
    println!("\n--- Time-Cached Approach ---");
    let start = Instant::now();
    let mut cached_count = 0;
    let mut time_cache: HashMap<DateTime<Utc>, spa::SpaTimeDependent> = HashMap::new();

    for &time in &times {
        let time_parts = time_cache
            .entry(time)
            .or_insert_with(|| spa::spa_time_dependent_parts(time, 69.0).unwrap());

        for &(lat, lon) in &coordinates {
            let _result = spa::spa_with_time_dependent_parts(
                time,
                lat,
                lon,
                0.0,
                69.0,
                Some(solar_positioning::RefractionCorrection::new(1013.25, 15.0).unwrap()),
                time_parts,
            )
            .unwrap();
            cached_count += 1;
        }
    }
    let cached_duration = start.elapsed();

    println!(
        "Completed {} calculations in {:.3}s",
        cached_count,
        cached_duration.as_secs_f64()
    );
    println!(
        "Rate: {:.0} calculations/second",
        cached_count as f64 / cached_duration.as_secs_f64()
    );
    println!("Cache entries used: {}", time_cache.len());

    // Results
    assert_eq!(naive_count, cached_count);
    assert_eq!(naive_count, total_calculations);

    let speedup = naive_duration.as_secs_f64() / cached_duration.as_secs_f64();

    println!("\n=== Results ===");
    println!("Naive: {:.3}s", naive_duration.as_secs_f64());
    println!("Cached: {:.3}s", cached_duration.as_secs_f64());
    println!("Speedup: {:.2}x", speedup);

    println!("\n=== Comparison to sunce ===");
    println!("sunce typically takes ~19s for this workload");
    println!(
        "Library naive: {:.1}s ({:.1}x faster than sunce)",
        naive_duration.as_secs_f64(),
        19.0 / naive_duration.as_secs_f64()
    );
    println!(
        "Library cached: {:.1}s ({:.1}x faster than sunce)",
        cached_duration.as_secs_f64(),
        19.0 / cached_duration.as_secs_f64()
    );
}
