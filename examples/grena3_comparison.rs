//! Compare SPA and Grena3 algorithms for speed vs accuracy trade-offs.

use chrono::{DateTime, Utc};
use solar_positioning::{grena3, spa, time::DeltaT};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>()?;
    let latitude = 37.7749; // San Francisco
    let longitude = -122.4194;
    let delta_t = DeltaT::estimate_from_date(2023, 6)?;

    println!("Algorithm comparison for San Francisco on June 21, 2023:");
    println!();

    // SPA calculation
    let start = Instant::now();
    let spa_position =
        spa::solar_position(datetime, latitude, longitude, 0.0, delta_t, 1013.25, 15.0)?;
    let spa_duration = start.elapsed();

    println!("SPA (High Accuracy):");
    println!("  Azimuth: {:.6}°", spa_position.azimuth());
    println!("  Elevation: {:.6}°", spa_position.elevation_angle());
    println!("  Time: {:?}", spa_duration);
    println!();

    // Grena3 calculation
    let start = Instant::now();
    let grena3_position = grena3::solar_position(datetime, latitude, longitude, delta_t)?;
    let grena3_duration = start.elapsed();

    println!("Grena3 (Fast):");
    println!("  Azimuth: {:.6}°", grena3_position.azimuth());
    println!("  Elevation: {:.6}°", grena3_position.elevation_angle());
    println!("  Time: {:?}", grena3_duration);
    println!();

    // Calculate differences
    let azimuth_diff = (spa_position.azimuth() - grena3_position.azimuth()).abs();
    let elevation_diff = (spa_position.elevation_angle() - grena3_position.elevation_angle()).abs();

    println!("Differences:");
    println!("  Azimuth: {:.6}°", azimuth_diff);
    println!("  Elevation: {:.6}°", elevation_diff);
    println!();

    if spa_duration.as_nanos() > 0 && grena3_duration.as_nanos() > 0 {
        let speedup = spa_duration.as_nanos() as f64 / grena3_duration.as_nanos() as f64;
        println!("Grena3 is {:.1}x faster than SPA", speedup);
    }

    println!();
    println!("Note: Grena3 is designed for years 2010-2110 with ~0.01° accuracy.");
    println!("SPA works for years -2000 to 6000 with ~0.0003° accuracy.");

    Ok(())
}
