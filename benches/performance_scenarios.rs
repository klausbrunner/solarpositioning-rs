use chrono::{DateTime, Duration, Utc};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use solar_positioning::{grena3, spa};
use std::hint::black_box;

/// Benchmark configuration for different usage patterns
#[derive(Clone)]
#[allow(dead_code)]
struct BenchmarkScenario {
    name: &'static str,
    description: &'static str,
}

#[allow(dead_code)]
const SCENARIOS: &[BenchmarkScenario] = &[
    BenchmarkScenario {
        name: "single_calculation",
        description: "Single solar position calculation (baseline)",
    },
    BenchmarkScenario {
        name: "time_series_fixed_location",
        description: "Time series at fixed location (weather station pattern)",
    },
    BenchmarkScenario {
        name: "coordinate_sweep_fixed_time",
        description: "Geographic grid at fixed time (solar resource mapping)",
    },
    BenchmarkScenario {
        name: "mixed_coordinates_and_times",
        description: "Combined coordinate and time variations",
    },
    BenchmarkScenario {
        name: "random_access_pattern",
        description: "Random locations and times (worst case for caching)",
    },
];

fn benchmark_single_calculation(c: &mut Criterion) {
    let mut group = c.benchmark_group("single");
    group.sample_size(10);
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.measurement_time(std::time::Duration::from_secs(3));

    let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let lat = 37.7749;
    let lon = -122.4194;

    group.bench_function("spa", |b| {
        b.iter(|| {
            spa::solar_position(
                black_box(datetime),
                black_box(lat),
                black_box(lon),
                black_box(0.0),
                black_box(69.0),
                black_box(Some(
                    solar_positioning::RefractionCorrection::new(1013.25, 15.0).unwrap(),
                )),
            )
            .unwrap()
        })
    });

    group.bench_function("grena3", |b| {
        b.iter(|| {
            grena3::solar_position(
                black_box(datetime),
                black_box(lat),
                black_box(lon),
                black_box(69.0),
                black_box(None),
            )
            .unwrap()
        })
    });

    group.finish();
}

fn benchmark_time_series_fixed_location(c: &mut Criterion) {
    let mut group = c.benchmark_group("time_series_fixed_location");
    group.sample_size(10);
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.measurement_time(std::time::Duration::from_secs(5));

    let base_datetime = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let lat = 37.7749;
    let lon = -122.4194;

    for &count in &[100_000, 1_000_000] {
        // 100K and 1M calculations
        group.throughput(Throughput::Elements(count));

        let datetimes: Vec<DateTime<Utc>> = (0..count)
            .map(|i| base_datetime + Duration::hours(i as i64))
            .collect();

        group.bench_with_input(BenchmarkId::new("spa", count), &count, |b, _| {
            b.iter(|| {
                for &dt in &datetimes {
                    let _result = spa::solar_position(
                        black_box(dt),
                        black_box(lat),
                        black_box(lon),
                        black_box(0.0),
                        black_box(69.0),
                        black_box(Some(
                            solar_positioning::RefractionCorrection::new(1013.25, 15.0).unwrap(),
                        )),
                    )
                    .unwrap();
                }
            })
        });

        group.bench_with_input(BenchmarkId::new("grena3", count), &count, |b, _| {
            b.iter(|| {
                for &dt in &datetimes {
                    let _result = grena3::solar_position(
                        black_box(dt),
                        black_box(lat),
                        black_box(lon),
                        black_box(69.0),
                        black_box(None),
                    )
                    .unwrap();
                }
            })
        });
    }

    group.finish();
}

fn benchmark_coordinate_sweep_fixed_time(c: &mut Criterion) {
    let mut group = c.benchmark_group("coordinate_sweep_fixed_time");
    group.sample_size(10);
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.measurement_time(std::time::Duration::from_secs(5));

    let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();

    for &grid_size in &[200, 1000] {
        // 200x200, 1000x1000 grids (40K, 1M calculations)
        let count = grid_size * grid_size;
        group.throughput(Throughput::Elements(count as u64));

        let coordinates: Vec<(f64, f64)> = (0..grid_size)
            .flat_map(|i| {
                (0..grid_size).map(move |j| {
                    let lat = -44.9 + (i as f64) * 89.8 / grid_size as f64; // -44.9° to +44.9° latitude (avoid ±90°)
                    let lon = -179.9 + (j as f64) * 359.8 / grid_size as f64; // Full longitude range
                    (lat, lon)
                })
            })
            .collect();

        group.bench_with_input(
            BenchmarkId::new("spa", format!("{}x{}", grid_size, grid_size)),
            &count,
            |b, _| {
                b.iter(|| {
                    for &(lat, lon) in &coordinates {
                        let _result = spa::solar_position(
                            black_box(datetime),
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
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("grena3", format!("{}x{}", grid_size, grid_size)),
            &count,
            |b, _| {
                b.iter(|| {
                    for &(lat, lon) in &coordinates {
                        let _result = grena3::solar_position(
                            black_box(datetime),
                            black_box(lat),
                            black_box(lon),
                            black_box(69.0),
                            black_box(None),
                        )
                        .unwrap();
                    }
                })
            },
        );
    }

    group.finish();
}

fn benchmark_mixed_coordinates_and_times(c: &mut Criterion) {
    let mut group = c.benchmark_group("mixed_coordinates_and_times");
    group.sample_size(10);
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.measurement_time(std::time::Duration::from_secs(5));

    let base_datetime = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();

    // 500 coords × 2000 times = 1M calculations (matches sunce pattern)
    let coords: Vec<(f64, f64)> = (0..500)
        .map(|i| {
            let lat = -44.9 + (i as f64) * 89.8 / 500.0;
            let lon = -179.9 + (i as f64) * 359.8 / 500.0;
            (lat, lon)
        })
        .collect();

    let times: Vec<DateTime<Utc>> = (0..2000)
        .map(|i| base_datetime + Duration::hours(i as i64))
        .collect();

    let total = (coords.len() * times.len()) as u64;
    group.throughput(Throughput::Elements(total));

    // Optimized: pre-compute time-dependent parts (sunce pattern)
    group.bench_function("spa_optimized", |b| {
        b.iter(|| {
            for &time in &times {
                let time_parts =
                    spa::spa_time_dependent_parts(black_box(time), black_box(69.0)).unwrap();
                for &(lat, lon) in &coords {
                    let _result = spa::spa_with_time_dependent_parts(
                        black_box(lat),
                        black_box(lon),
                        black_box(0.0),
                        black_box(Some(
                            solar_positioning::RefractionCorrection::new(1013.25, 15.0).unwrap(),
                        )),
                        &time_parts,
                    )
                    .unwrap();
                }
            }
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_single_calculation,
    benchmark_time_series_fixed_location,
    benchmark_coordinate_sweep_fixed_time,
    benchmark_mixed_coordinates_and_times
);

criterion_main!(benches);
