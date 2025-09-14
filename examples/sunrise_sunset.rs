//! Sunrise/sunset calculation example with different twilight types across diverse global locations.

use chrono::{DateTime, NaiveDate, Utc};
use solar_positioning::{Horizon, spa, time::DeltaT, types::SunriseResult};

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
    let date = NaiveDate::from_ymd_opt(2023, 12, 21)
        .unwrap()
        .and_hms_opt(0, 0, 0)
        .unwrap()
        .and_utc();

    let delta_t = DeltaT::estimate_from_date(2023, 12)?;

    for city in &cities {
        println!("=== {} ===", city.name);
        println!(
            "Coordinates: {:.2}°N, {:.2}°E",
            city.latitude, city.longitude
        );
        println!("Date: December 21, 2023 (Winter Solstice)");
        println!();

        calculate_and_print_times(date, city.latitude, city.longitude, delta_t)?;
        println!();
    }

    Ok(())
}

fn calculate_and_print_times(
    date: DateTime<Utc>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    // Calculate all twilight types
    let horizons = [
        ("Sunrise/Sunset", Horizon::SunriseSunset),
        ("Civil Twilight", Horizon::CivilTwilight),
        ("Nautical Twilight", Horizon::NauticalTwilight),
        ("Astronomical Twilight", Horizon::AstronomicalTwilight),
    ];

    for (name, horizon) in &horizons {
        let result = spa::sunrise_sunset_for_horizon(date, latitude, longitude, delta_t, *horizon)?;
        print_sunrise_result(name, &result);
    }

    Ok(())
}

fn print_sunrise_result(label: &str, result: &SunriseResult<DateTime<Utc>>) {
    println!("{}:", label);
    match result {
        SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        } => {
            println!("  Begin: {}", sunrise.format("%H:%M:%S UTC"));
            println!("  Transit: {}", transit.format("%H:%M:%S UTC"));
            println!("  End: {}", sunset.format("%H:%M:%S UTC"));
        }
        SunriseResult::AllDay { transit } => {
            println!("  All day above horizon");
            println!("  Transit: {}", transit.format("%H:%M:%S UTC"));
        }
        SunriseResult::AllNight { transit } => {
            println!("  All night below horizon");
            println!("  Transit: {}", transit.format("%H:%M:%S UTC"));
        }
    }
    println!();
}
