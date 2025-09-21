//! Test combined coordinate+time optimization patterns.

#![cfg(feature = "unstable")]

use chrono::{DateTime, Duration, Utc};
use solar_positioning::spa;
use std::collections::HashMap;
use std::time::Instant;

#[test]
fn test_combined_pattern_realistic_scenario() {
    // Realistic scenario: weather model covering small region over 24 hours
    let base_datetime = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();

    // 10x10 grid over southern Germany (realistic for regional weather models)
    let coordinates: Vec<(f64, f64)> = (0..10)
        .flat_map(|i| {
            (0..10).map(move |j| {
                let lat = 47.0 + (i as f64) * 0.5; // 47°-51.5° N (southern Germany)
                let lon = 7.0 + (j as f64) * 0.5; // 7°-11.5° E
                (lat, lon)
            })
        })
        .collect();

    // 24 hourly time steps
    let times: Vec<DateTime<Utc>> = (0..24)
        .map(|h| base_datetime + Duration::hours(h))
        .collect();

    println!(
        "Testing {}×{} grid over {} hours = {} calculations",
        10,
        10,
        times.len(),
        coordinates.len() * times.len()
    );

    // Strategy 1: Naive approach (full SPA for each coordinate+time pair)
    let start = Instant::now();
    let mut naive_results = Vec::new();

    for &time in &times {
        for &(lat, lon) in &coordinates {
            let result = spa::solar_position(time, lat, lon, 100.0, 69.0, 1013.25, 15.0).unwrap();
            naive_results.push(result);
        }
    }
    let naive_duration = start.elapsed();

    // Strategy 2: Time-based caching (cache time-dependent parts per time)
    let start = Instant::now();
    let mut cached_results = Vec::new();
    let mut time_cache: HashMap<DateTime<Utc>, spa::SpaTimeDependent> = HashMap::new();

    for &time in &times {
        // Get or calculate time-dependent parts for this time
        let time_parts = time_cache
            .entry(time)
            .or_insert_with(|| spa::spa_time_dependent_parts(time, 69.0).unwrap());

        // Fast calculation for each coordinate using cached time parts
        for &(lat, lon) in &coordinates {
            let result = spa::spa_with_time_dependent_parts(
                time, lat, lon, 100.0, 69.0, 1013.25, 15.0, time_parts,
            )
            .unwrap();
            cached_results.push(result);
        }
    }
    let cached_duration = start.elapsed();

    // Verify results are identical
    assert_eq!(naive_results.len(), cached_results.len());
    for (naive, cached) in naive_results.iter().zip(cached_results.iter()) {
        assert!((naive.azimuth() - cached.azimuth()).abs() < 1e-10);
        assert!((naive.zenith_angle() - cached.zenith_angle()).abs() < 1e-10);
    }

    let speedup = naive_duration.as_secs_f64() / cached_duration.as_secs_f64();

    println!("Naive approach: {:.3}s", naive_duration.as_secs_f64());
    println!(
        "Time-cached approach: {:.3}s",
        cached_duration.as_secs_f64()
    );
    println!("Speedup: {:.2}x", speedup);
    println!("Time cache entries: {}", time_cache.len());

    // Validate that we get meaningful speedup for this combined pattern
    assert!(
        speedup > 1.5,
        "Expected speedup > 1.5x, got {:.2}x",
        speedup
    );
}

#[test]
fn test_combined_pattern_varying_time_density() {
    // Test how speedup changes with different time vs coordinate densities
    let base_datetime = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();

    let test_scenarios = vec![
        (5, 2),  // 5 coords, 2 times (high coord density)
        (5, 10), // 5 coords, 10 times (balanced)
        (2, 10), // 2 coords, 10 times (high time density)
        (10, 1), // 10 coords, 1 time (pure coordinate sweep)
        (1, 10), // 1 coord, 10 times (pure time series)
    ];

    for (coord_count, time_count) in test_scenarios {
        // Generate coordinates
        let coordinates: Vec<(f64, f64)> = (0..coord_count)
            .map(|i| (50.0 + i as f64, 10.0 + i as f64))
            .collect();

        // Generate times
        let times: Vec<DateTime<Utc>> = (0..time_count)
            .map(|h| base_datetime + Duration::hours(h))
            .collect();

        // Time-cached approach
        let start = Instant::now();
        let mut time_cache: HashMap<DateTime<Utc>, spa::SpaTimeDependent> = HashMap::new();
        let mut result_count = 0;

        for &time in &times {
            let time_parts = time_cache
                .entry(time)
                .or_insert_with(|| spa::spa_time_dependent_parts(time, 69.0).unwrap());

            for &(lat, lon) in &coordinates {
                let _result = spa::spa_with_time_dependent_parts(
                    time, lat, lon, 0.0, 69.0, 1013.25, 15.0, time_parts,
                )
                .unwrap();
                result_count += 1;
            }
        }
        let cached_duration = start.elapsed();

        println!(
            "{}C×{}T ({} calcs): {:.3}s, {} cache entries",
            coord_count,
            time_count,
            result_count,
            cached_duration.as_secs_f64(),
            time_cache.len()
        );

        // Verify we calculated expected number of results
        assert_eq!(result_count, coord_count * time_count);
        assert_eq!(time_cache.len(), time_count as usize);
    }
}

#[test]
fn test_combined_vs_pure_patterns() {
    // Compare combined patterns against pure coordinate sweeps and pure time series
    let base_datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();

    // Pure coordinate sweep (100 coordinates, 1 time)
    let coordinates: Vec<(f64, f64)> = (0..100)
        .map(|i| (45.0 + (i as f64) * 0.1, 8.0 + (i as f64) * 0.1))
        .collect();

    println!("\n=== Pure Coordinate Sweep (100 coords, 1 time) ===");
    let start = Instant::now();
    let time_parts = spa::spa_time_dependent_parts(base_datetime, 69.0).unwrap();
    for &(lat, lon) in &coordinates {
        let _result = spa::spa_with_time_dependent_parts(
            base_datetime,
            lat,
            lon,
            0.0,
            69.0,
            1013.25,
            15.0,
            &time_parts,
        )
        .unwrap();
    }
    let coord_sweep_duration = start.elapsed();
    println!("Duration: {:.3}s", coord_sweep_duration.as_secs_f64());

    // Combined pattern (10 coordinates, 10 times = 100 calculations)
    let coords_small: Vec<(f64, f64)> = coordinates[0..10].to_vec();
    let times: Vec<DateTime<Utc>> = (0..10)
        .map(|h| base_datetime + Duration::hours(h))
        .collect();

    println!("\n=== Combined Pattern (10 coords, 10 times) ===");
    let start = Instant::now();
    let mut time_cache: HashMap<DateTime<Utc>, spa::SpaTimeDependent> = HashMap::new();
    for &time in &times {
        let time_parts = time_cache
            .entry(time)
            .or_insert_with(|| spa::spa_time_dependent_parts(time, 69.0).unwrap());
        for &(lat, lon) in &coords_small {
            let _result = spa::spa_with_time_dependent_parts(
                time, lat, lon, 0.0, 69.0, 1013.25, 15.0, time_parts,
            )
            .unwrap();
        }
    }
    let combined_duration = start.elapsed();
    println!("Duration: {:.3}s", combined_duration.as_secs_f64());
    println!("Time cache entries: {}", time_cache.len());

    // Analysis
    println!("\n=== Analysis ===");
    println!("Both scenarios: 100 total calculations");
    println!("Pure coordinate sweep: leverages maximum time sharing");
    println!("Combined pattern: partial time sharing (10 time points)");

    let efficiency_ratio = coord_sweep_duration.as_secs_f64() / combined_duration.as_secs_f64();
    println!(
        "Efficiency ratio (coord_sweep/combined): {:.2}",
        efficiency_ratio
    );

    // Combined should be slower than pure coordinate sweep due to less sharing
    // but still faster than naive approach
}
