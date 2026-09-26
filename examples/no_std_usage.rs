//! Numeric positions without chrono. The library calls also work without std or allocation.

use solar_positioning::{Location, RefractionCorrection, SolarPositions, time::JulianDate};

fn main() -> solar_positioning::Result<()> {
    let time = JulianDate::from_utc(2024, 6, 21, 12, 0, 0.0, 69.0)?;
    let atmosphere = Some(RefractionCorrection::standard());
    let locations = [
        (
            "Vienna",
            Location {
                latitude: 48.21,
                longitude: 16.37,
            },
        ),
        (
            "San Francisco",
            Location {
                latitude: 37.7749,
                longitude: -122.4194,
            },
        ),
        (
            "Sydney",
            Location {
                latitude: -33.8688,
                longitude: 151.2093,
            },
        ),
    ];
    for (name, positions) in [
        ("SPA", SolarPositions::new()),
        ("Grena3", SolarPositions::grena3()),
    ] {
        let prepared = positions.for_time_from_julian(time);
        for (city, location) in locations {
            let position = prepared.at(location, 0.0, atmosphere)?;
            println!(
                "{name}, {city}: azimuth {:.3}°, elevation {:.3}°",
                position.azimuth(),
                position.elevation_angle()
            );
        }
    }
    Ok(())
}
