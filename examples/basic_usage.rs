//! Basic solar position calculation example.

use chrono::{DateTime, Utc};
use solar_positioning::{spa, time::DeltaT};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Calculate solar position for San Francisco on summer solstice
    let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>()?;
    let latitude = 37.7749; // San Francisco
    let longitude = -122.4194;

    // Use estimated DeltaT for 2023
    let delta_t = DeltaT::estimate_from_date(2023, 6)?;

    // Calculate position with atmospheric refraction correction
    let position = spa::solar_position(
        datetime, latitude, longitude, 0.0,     // elevation (m)
        delta_t, // deltaT (s)
        1013.25, // pressure (hPa)
        15.0,    // temperature (°C)
    )?;

    println!("Solar position for San Francisco on June 21, 2023 at noon UTC:");
    println!("  Azimuth: {:.3}°", position.azimuth());
    println!("  Elevation: {:.3}°", position.elevation_angle());
    println!("  Zenith angle: {:.3}°", position.zenith_angle());

    if position.is_sun_up() {
        println!("  Sun is above the horizon");
    } else {
        println!("  Sun is below the horizon");
    }

    Ok(())
}
