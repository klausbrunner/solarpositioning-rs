//! Compare positions from SPA and Grena3 with the same inputs.

use chrono::{DateTime, Utc};
use solar_positioning::{Location, SolarPositions, delta_t};

fn main() -> solar_positioning::Result<()> {
    let time = "2023-06-21T19:00:00Z".parse::<DateTime<Utc>>().unwrap();
    let location = Location {
        latitude: 37.7749,
        longitude: -122.4194,
    };
    let delta_t = delta_t::estimate_from_date_like(time.date_naive())?;
    for (name, positions) in [
        ("SPA", SolarPositions::new()),
        ("Grena3", SolarPositions::grena3()),
    ] {
        let position = positions.at(&time, location, 0.0, delta_t, None)?;
        println!(
            "{name}: azimuth {:.6}°, elevation {:.6}°",
            position.azimuth(),
            position.elevation_angle()
        );
    }
    Ok(())
}
