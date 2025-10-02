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
    let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let lat = 37.7749;
    let lon = -122.4194;

    c.bench_function("spa_single", |b| {
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

    c.bench_function("grena3_single", |b| {
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
}

fn benchmark_time_series_fixed_location(c: &mut Criterion) {
    let mut group = c.benchmark_group("time_series_fixed_location");

    let base_datetime = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let lat = 37.7749;
    let lon = -122.4194;

    for &count in &[1000, 5000, 25000] {
        // Sized for ~5-15 second runtimes
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

    let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();

    for &grid_size in &[30, 70, 150] {
        // 30x30, 70x70, 150x150 grids (~1K, 5K, 22K calculations)
        let count = grid_size * grid_size;
        group.throughput(Throughput::Elements(count as u64));

        let coordinates: Vec<(f64, f64)> = (0..grid_size)
            .flat_map(|i| {
                (0..grid_size).map(move |j| {
                    let lat = 30.0 + (i as f64) * 0.1; // 30° to 40° latitude
                    let lon = -120.0 + (j as f64) * 0.1; // -120° to -110° longitude
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

    let base_datetime = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();

    for &(coord_count, time_count) in &[(20, 50), (50, 100), (100, 250)] {
        // ~1K, 5K, 25K calculations
        let total_count = coord_count * time_count;
        group.throughput(Throughput::Elements(total_count as u64));

        let test_data: Vec<(f64, f64, DateTime<Utc>)> = (0..coord_count)
            .flat_map(|i| {
                (0..time_count).map(move |j| {
                    let lat = 30.0 + (i as f64) * 0.2;
                    let lon = -120.0 + (i as f64) * 0.2;
                    let dt = base_datetime + Duration::hours(j as i64);
                    (lat, lon, dt)
                })
            })
            .collect();

        group.bench_with_input(
            BenchmarkId::new("spa", format!("{}coords_{}times", coord_count, time_count)),
            &total_count,
            |b, _| {
                b.iter(|| {
                    for &(lat, lon, dt) in &test_data {
                        let _result = spa::solar_position(
                            black_box(dt),
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
    }

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
