//! Solar position with a timezone-aware timestamp.

use chrono::{DateTime, FixedOffset, Utc};
use solar_positioning::{Location, RefractionCorrection, SolarPositions, delta_t};

fn main() -> solar_positioning::Result<()> {
    let time = "2023-06-21T12:00:00-07:00"
        .parse::<DateTime<FixedOffset>>()
        .unwrap();
    let location = Location {
        latitude: 37.7749,
        longitude: -122.4194,
    }; // San Francisco
    let delta_t = delta_t::estimate_from_date_like(time.date_naive())?;
    let atmosphere = Some(RefractionCorrection::standard());
    let positions = SolarPositions::new();
    let position = positions.at(&time, location, 0.0, delta_t, atmosphere)?;

    println!("San Francisco at {time}:");
    println!("Azimuth: {:.3}°", position.azimuth());
    println!("Elevation: {:.3}°", position.elevation_angle());
    assert_eq!(
        position,
        positions.at(
            &time.with_timezone(&Utc),
            location,
            0.0,
            delta_t,
            atmosphere
        )?
    );
    Ok(())
}
