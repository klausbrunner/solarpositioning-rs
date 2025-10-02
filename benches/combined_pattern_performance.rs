use chrono::{DateTime, Duration, Utc};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use solar_positioning::spa;
use std::collections::HashMap;
use std::hint::black_box;

fn benchmark_combined_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("combined_patterns");

    let base_datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();

    // Test different combinations of coordinates and times
    let test_cases = vec![
        (10, 10), // 10×10 = 100 calculations (balanced)
        (50, 5),  // 50×5 = 250 calculations (coord-heavy)
        (5, 50),  // 5×50 = 250 calculations (time-heavy)
        (100, 1), // 100×1 = 100 calculations (pure coord sweep)
        (1, 100), // 1×100 = 100 calculations (pure time series)
        (25, 20), // 25×20 = 500 calculations (realistic weather model)
    ];

    for (coord_count, time_count) in test_cases {
        let total_calculations = coord_count * time_count;
        group.throughput(Throughput::Elements(total_calculations as u64));

        // Generate test data
        let coordinates: Vec<(f64, f64)> = (0..coord_count)
            .map(|i| {
                let lat = 45.0 + (i as f64) * 0.1;
                let lon = 8.0 + (i as f64) * 0.1;
                (lat, lon)
            })
            .collect();

        let times: Vec<DateTime<Utc>> = (0..time_count)
            .map(|h| base_datetime + Duration::hours(h))
            .collect();

        // Naive approach
        group.bench_with_input(
            BenchmarkId::new("naive", format!("{}C×{}T", coord_count, time_count)),
            &total_calculations,
            |b, _| {
                b.iter(|| {
                    for &time in &times {
                        for &(lat, lon) in &coordinates {
                            let _result = spa::solar_position(
                                black_box(time),
                                black_box(lat),
                                black_box(lon),
                                black_box(0.0),
                                black_box(69.0),
                                black_box(Some(
                                    solar_positioning::RefractionCorrection::new(1013.25, 15.0)
                                        .unwrap(),
                                )),
                            )
                            .unwrap();
                        }
                    }
                })
            },
        );

        // Time-cached approach
        group.bench_with_input(
            BenchmarkId::new("time_cached", format!("{}C×{}T", coord_count, time_count)),
            &total_calculations,
            |b, _| {
                b.iter(|| {
                    let mut time_cache: HashMap<DateTime<Utc>, spa::SpaTimeDependent> =
                        HashMap::new();
                    for &time in &times {
                        let time_parts = time_cache.entry(time).or_insert_with(|| {
                            spa::spa_time_dependent_parts(black_box(time), black_box(69.0)).unwrap()
                        });

                        for &(lat, lon) in &coordinates {
                            let _result = spa::spa_with_time_dependent_parts(
                                black_box(lat),
                                black_box(lon),
                                black_box(0.0),
                                black_box(Some(
                                    solar_positioning::RefractionCorrection::new(1013.25, 15.0)
                                        .unwrap(),
                                )),
                                black_box(time_parts),
                            )
                            .unwrap();
                        }
                    }
                })
            },
        );
    }

    group.finish();
}

fn benchmark_optimization_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("optimization_scaling");

    let base_datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();

    // Test how optimization effectiveness changes with time/coordinate ratio
    for time_count in [1, 5, 10, 25, 50] {
        let coord_count = 100; // Fixed coordinate count
        let total = coord_count * time_count;

        group.throughput(Throughput::Elements(total as u64));

        let coordinates: Vec<(f64, f64)> = (0..coord_count)
            .map(|i| (50.0 + (i as f64) * 0.01, 10.0 + (i as f64) * 0.01))
            .collect();

        let times: Vec<DateTime<Utc>> = (0..time_count)
            .map(|h| base_datetime + Duration::hours(h))
            .collect();

        group.bench_with_input(
            BenchmarkId::new("optimized", format!("{}T", time_count)),
            &total,
            |b, _| {
                b.iter(|| {
                    let mut time_cache: HashMap<DateTime<Utc>, spa::SpaTimeDependent> =
                        HashMap::new();
                    for &time in &times {
                        let time_parts = time_cache.entry(time).or_insert_with(|| {
                            spa::spa_time_dependent_parts(black_box(time), black_box(69.0)).unwrap()
                        });

                        for &(lat, lon) in &coordinates {
                            let _result = spa::spa_with_time_dependent_parts(
                                black_box(lat),
                                black_box(lon),
                                black_box(0.0),
                                black_box(Some(
                                    solar_positioning::RefractionCorrection::new(1013.25, 15.0)
                                        .unwrap(),
                                )),
                                black_box(time_parts),
                            )
                            .unwrap();
                        }
                    }
                })
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_combined_patterns,
    benchmark_optimization_scaling
);

criterion_main!(benches);
