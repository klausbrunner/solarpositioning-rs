//! Sunrise/sunset calculation example with different twilight types across diverse global locations.

use chrono::{DateTime, TimeZone, Utc};
use solar_positioning::{spa, time::DeltaT, types::SunriseResult, Horizon};

#[derive(Debug)]
struct City {
    name: &'static str,
    latitude: f64,
    longitude: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Test cities from around the world (diverse latitudes and longitudes)
    let cities = [
        City {
            name: "Longyearbyen, Norway (Arctic)",
            latitude: 78.22,
            longitude: 15.65,
        },
        City {
            name: "Anchorage, Alaska",
            latitude: 61.216667,
            longitude: -149.866667,
        },
        City {
            name: "Auckland, New Zealand",
            latitude: -36.840556,
            longitude: 174.74,
        },
        City {
            name: "Singapore",
            latitude: 1.283333,
            longitude: 103.833333,
        },
        City {
            name: "Brasília, Brazil",
            latitude: -15.8,
            longitude: -47.85,
        },
    ];

    // Calculate for winter solstice - this shows the most extreme variations
    let date_utc = Utc.with_ymd_and_hms(2023, 12, 21, 0, 0, 0).unwrap();
    let delta_t = DeltaT::estimate_from_date(2023, 12)?;

    for city in &cities {
        println!("=== {} ===", city.name);
        println!(
            "Coordinates: {:.2}°N, {:.2}°E",
            city.latitude, city.longitude
        );
        println!("Date: December 21, 2023 (Winter Solstice)");
        println!();

        calculate_and_print_times(date_utc, city.latitude, city.longitude, delta_t)?;
        println!();
    }

    Ok(())
}

fn calculate_and_print_times<Tz: TimeZone>(
    date: DateTime<Tz>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
) -> Result<(), Box<dyn std::error::Error>>
where
    <Tz as TimeZone>::Offset: std::fmt::Display,
{
    // Calculate all twilight types
    let horizons = [
        ("Sunrise/Sunset", Horizon::SunriseSunset),
        ("Civil Twilight", Horizon::CivilTwilight),
        ("Nautical Twilight", Horizon::NauticalTwilight),
        ("Astronomical Twilight", Horizon::AstronomicalTwilight),
    ];

    for (name, horizon) in &horizons {
        let result =
            spa::sunrise_sunset_for_horizon(date.clone(), latitude, longitude, delta_t, *horizon)?;
        print_sunrise_result(name, &result);
    }

    Ok(())
}

fn print_sunrise_result<Tz: TimeZone>(label: &str, result: &SunriseResult<DateTime<Tz>>)
where
    <Tz as TimeZone>::Offset: std::fmt::Display,
{
    println!("{}:", label);
    match result {
        SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        } => {
            println!("  Begin: {}", sunrise.format("%H:%M:%S %z"));
            println!("  Transit: {}", transit.format("%H:%M:%S %z"));
            println!("  End: {}", sunset.format("%H:%M:%S %z"));
        }
        SunriseResult::AllDay { transit } => {
            println!("  All day above horizon");
            println!("  Transit: {}", transit.format("%H:%M:%S %z"));
        }
        SunriseResult::AllNight { transit } => {
            println!("  All night below horizon");
            println!("  Transit: {}", transit.format("%H:%M:%S %z"));
        }
    }
    println!();
}
