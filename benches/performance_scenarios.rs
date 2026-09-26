use chrono::{DateTime, Duration, Utc};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use solar_positioning::{Location, RefractionCorrection, SolarPositions, time::JulianDate};
use std::hint::black_box;

const MODELS: [(&str, SolarPositions); 2] = [
    ("spa", SolarPositions::new()),
    ("grena3", SolarPositions::grena3()),
];
const LOCATION: Location = Location {
    latitude: 37.7749,
    longitude: -122.4194,
};
const REFRACTION: Option<RefractionCorrection> = Some(RefractionCorrection::standard());

fn benchmark_single_calculation(c: &mut Criterion) {
    let mut group = c.benchmark_group("single");
    let time = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
    for (name, positions) in MODELS {
        group.bench_function(name, |b| {
            b.iter(|| {
                positions
                    .at(
                        black_box(&time),
                        black_box(LOCATION),
                        black_box(0.0),
                        black_box(69.0),
                        black_box(REFRACTION),
                    )
                    .unwrap()
            })
        });
    }
    group.finish();
}

fn benchmark_time_series_fixed_location(c: &mut Criterion) {
    let mut group = c.benchmark_group("time_series_fixed_location");
    group.sample_size(10);
    let start = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    for count in [100_000_u32, 1_000_000] {
        group.throughput(Throughput::Elements(u64::from(count)));
        // Half-hour steps keep even the longest series within both models' validity ranges.
        let times: Vec<_> = (0..count)
            .map(|i| start + Duration::minutes(i64::from(i) * 30))
            .collect();
        for (name, positions) in MODELS {
            group.bench_with_input(BenchmarkId::new(name, count), &times, |b, times| {
                b.iter(|| {
                    for time in times {
                        black_box(
                            positions
                                .at(black_box(time), black_box(LOCATION), 0.0, 69.0, REFRACTION)
                                .unwrap(),
                        );
                    }
                })
            });
        }
    }
    group.finish();
}

fn benchmark_coordinate_sweep_fixed_time(c: &mut Criterion) {
    let mut group = c.benchmark_group("coordinate_sweep_fixed_time");
    group.sample_size(10);
    let time = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
    for grid_size in [200_u32, 1000] {
        let locations: Vec<_> = (0..grid_size)
            .flat_map(|i| {
                (0..grid_size).map(move |j| Location {
                    latitude: -44.9 + f64::from(i) * 89.8 / f64::from(grid_size),
                    longitude: -179.9 + f64::from(j) * 359.8 / f64::from(grid_size),
                })
            })
            .collect();
        group.throughput(Throughput::Elements(locations.len() as u64));
        for (name, positions) in MODELS {
            let prepared = positions.for_time(&time, 69.0).unwrap();
            group.bench_with_input(
                BenchmarkId::new(name, grid_size),
                &locations,
                |b, locations| {
                    b.iter(|| {
                        for location in locations {
                            black_box(
                                positions
                                    .at(
                                        black_box(&time),
                                        black_box(*location),
                                        0.0,
                                        69.0,
                                        REFRACTION,
                                    )
                                    .unwrap(),
                            );
                        }
                    })
                },
            );
            group.bench_with_input(
                BenchmarkId::new(format!("{name}_prepared"), grid_size),
                &locations,
                |b, locations| {
                    b.iter(|| {
                        for location in locations {
                            black_box(prepared.at(black_box(*location), 0.0, REFRACTION).unwrap());
                        }
                    })
                },
            );
        }
    }
    group.finish();
}

fn benchmark_solar_events_multiple(c: &mut Criterion) {
    let mut group = c.benchmark_group("solar_events_multiple");
    group.sample_size(20);
    group.warm_up_time(std::time::Duration::from_secs(1));
    group.measurement_time(std::time::Duration::from_secs(4));

    let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let lat = 37.7749;
    let lon = -122.4194;

    let horizon_sets: &[(&str, &[solar_positioning::Horizon])] = &[
        ("1", &[solar_positioning::Horizon::SunriseSunset]),
        (
            "3",
            &[
                solar_positioning::Horizon::SunriseSunset,
                solar_positioning::Horizon::CivilTwilight,
                solar_positioning::Horizon::NauticalTwilight,
            ],
        ),
        (
            "5",
            &[
                solar_positioning::Horizon::SunriseSunset,
                solar_positioning::Horizon::CivilTwilight,
                solar_positioning::Horizon::NauticalTwilight,
                solar_positioning::Horizon::AstronomicalTwilight,
                solar_positioning::Horizon::Custom(-4.0),
            ],
        ),
    ];

    for &(name, horizons) in horizon_sets {
        group.throughput(Throughput::Elements(horizons.len() as u64));
        group.bench_with_input(BenchmarkId::new("spa", name), horizons, |b, horizons| {
            b.iter(|| {
                solar_positioning::SolarEvents::new()
                    .for_date_multiple(
                        black_box(datetime.date_naive()),
                        &Utc,
                        Location {
                            latitude: black_box(lat),
                            longitude: black_box(lon),
                        },
                        black_box(69.0),
                        horizons.iter().copied(),
                    )
                    .unwrap()
            })
        });
    }

    group.finish();
}

fn benchmark_preparation(c: &mut Criterion) {
    let mut group = c.benchmark_group("preparation");
    let time = JulianDate::from_utc(2023, 6, 21, 12, 0, 0.0, 69.0).unwrap();
    for (name, positions) in MODELS {
        group.bench_function(name, |b| {
            b.iter(|| positions.for_time_from_julian(black_box(time)))
        });
    }
    group.finish();
}

fn benchmark_mixed_coordinates_and_times(c: &mut Criterion) {
    let mut group = c.benchmark_group("mixed_coordinates_and_times");
    group.sample_size(10);
    let start = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();
    // 500 locations × 2000 timestamps, preparing once per timestamp.
    let locations: Vec<_> = (0..500)
        .map(|i| Location {
            latitude: -44.9 + f64::from(i) * 89.8 / 500.0,
            longitude: -179.9 + f64::from(i) * 359.8 / 500.0,
        })
        .collect();
    let times: Vec<_> = (0..2000).map(|i| start + Duration::hours(i)).collect();
    group.throughput(Throughput::Elements((locations.len() * times.len()) as u64));
    for (name, positions) in MODELS {
        group.bench_function(name, |b| {
            b.iter(|| {
                for time in &times {
                    let prepared = positions.for_time(black_box(time), 69.0).unwrap();
                    for location in &locations {
                        black_box(prepared.at(black_box(*location), 0.0, REFRACTION).unwrap());
                    }
                }
            })
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    benchmark_single_calculation,
    benchmark_time_series_fixed_location,
    benchmark_coordinate_sweep_fixed_time,
    benchmark_solar_events_multiple,
    benchmark_preparation,
    benchmark_mixed_coordinates_and_times
);
criterion_main!(benches);
