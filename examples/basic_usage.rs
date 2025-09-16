//! Basic solar position calculation example.

use chrono::{DateTime, FixedOffset, TimeZone, Utc};
use solar_positioning::{spa, time::DeltaT};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Example 1: Calculate solar position using FixedOffset timezone
    let datetime_fixed = "2023-06-21T12:00:00-07:00".parse::<DateTime<FixedOffset>>()?;

    // Example 2: Same time using UTC
    let datetime_utc = Utc.with_ymd_and_hms(2023, 6, 21, 19, 0, 0).unwrap(); // 19:00 UTC = 12:00 PDT
    let latitude = 37.7749; // San Francisco
    let longitude = -122.4194;

    // Use estimated DeltaT for 2023
    let delta_t = DeltaT::estimate_from_date(2023, 6)?;

    // Calculate position with atmospheric refraction correction using FixedOffset
    let position_fixed = spa::solar_position(
        datetime_fixed,
        latitude,
        longitude,
        0.0,     // elevation (m)
        delta_t, // deltaT (s)
        1013.25, // pressure (hPa)
        15.0,    // temperature (°C)
    )?;

    // Same calculation using UTC timezone
    let position_utc = spa::solar_position(
        datetime_utc,
        latitude,
        longitude,
        0.0,     // elevation (m)
        delta_t, // deltaT (s)
        1013.25, // pressure (hPa)
        15.0,    // temperature (°C)
    )?;

    println!("Solar position for San Francisco on June 21, 2023 at noon Pacific Time:");
    println!("Using FixedOffset timezone:");
    println!("  Azimuth: {:.3}°", position_fixed.azimuth());
    println!("  Elevation: {:.3}°", position_fixed.elevation_angle());
    println!("  Zenith angle: {:.3}°", position_fixed.zenith_angle());

    println!("\nUsing UTC timezone (same moment):");
    println!("  Azimuth: {:.3}°", position_utc.azimuth());
    println!("  Elevation: {:.3}°", position_utc.elevation_angle());
    println!("  Zenith angle: {:.3}°", position_utc.zenith_angle());

    println!(
        "\nBoth calculations produce identical results: {}",
        position_fixed.azimuth() == position_utc.azimuth()
            && position_fixed.elevation_angle() == position_utc.elevation_angle()
    );

    if position_fixed.is_sun_up() {
        println!("  Sun is above the horizon");
    } else {
        println!("  Sun is below the horizon");
    }

    Ok(())
}
