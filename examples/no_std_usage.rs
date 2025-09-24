//! Example demonstrating usage without std/chrono dependencies.
//!
//! This shows how to use the library in no_std environments where
//! users handle their own time conversions.

use solar_positioning::{RefractionCorrection, grena3, spa, time::JulianDate};

fn main() {
    // Example: Calculate solar position for 2024-06-21 12:00:00 UTC
    // Vienna: 48.21°N, 16.37°E

    println!("Solar positioning without std/chrono dependencies\n");

    // Create Julian date from UTC components
    let jd = JulianDate::from_utc(
        2024, // year
        6,    // month
        21,   // day
        12,   // hour (UTC)
        0,    // minute
        0.0,  // second
        69.0, // deltaT
    )
    .expect("Valid date");

    println!("Julian Date: {:.6}", jd.julian_date());
    println!("Delta T: {:.1} seconds\n", jd.delta_t());

    // Calculate position using SPA (high accuracy)
    let spa_position = spa::solar_position_from_julian(
        jd,
        48.21, // latitude
        16.37, // longitude
        190.0, // elevation (meters)
        Some(RefractionCorrection::standard()),
    )
    .expect("Valid coordinates");

    println!("SPA Results:");
    println!("  Azimuth: {:.3}°", spa_position.azimuth());
    println!("  Elevation: {:.3}°", spa_position.elevation_angle());
    println!("  Zenith: {:.3}°\n", spa_position.zenith_angle());

    // Calculate using Grena3 (faster, less accurate)
    // First calculate the t parameter for Grena3
    let t = grena3::calc_t_from_components(2024, 6, 21, 12, 0, 0.0);

    let grena3_position = grena3::solar_position_from_t(
        t,
        48.21, // latitude
        16.37, // longitude
        69.0,  // deltaT
        Some(RefractionCorrection::standard()),
    )
    .expect("Valid coordinates");

    println!("Grena3 Results:");
    println!("  Azimuth: {:.3}°", grena3_position.azimuth());
    println!("  Elevation: {:.3}°", grena3_position.elevation_angle());
    println!("  Zenith: {:.3}°\n", grena3_position.zenith_angle());

    // Show difference between algorithms
    let azimuth_diff = (spa_position.azimuth() - grena3_position.azimuth()).abs();
    let elevation_diff = (spa_position.elevation_angle() - grena3_position.elevation_angle()).abs();

    println!("Algorithm Differences:");
    println!("  Azimuth difference: {:.4}°", azimuth_diff);
    println!("  Elevation difference: {:.4}°", elevation_diff);

    // Example: Using pre-computed time-dependent parts for coordinate sweeps
    println!("\nCoordinate sweep example (3 locations, same time):");

    let time_parts = spa::spa_time_dependent_from_julian(jd).expect("Valid Julian date");

    let locations = [
        ("Vienna", 48.21, 16.37),
        ("San Francisco", 37.7749, -122.4194),
        ("Sydney", -33.8688, 151.2093),
    ];

    for (name, lat, lon) in &locations {
        let position = spa::spa_with_time_dependent_parts(
            *lat,
            *lon,
            0.0, // sea level
            Some(RefractionCorrection::standard()),
            &time_parts,
        )
        .expect("Valid coordinates");

        println!(
            "  {} - Azimuth: {:.1}°, Elevation: {:.1}°",
            name,
            position.azimuth(),
            position.elevation_angle()
        );
    }
}
