//! Example demonstrating sunrise/sunset calculation without the chrono library.
//!
//! This example shows how to use the numeric API when you don't want to depend on chrono.

use solar_positioning::{spa, Horizon};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Calculate sunrise/sunset for San Francisco on June 21, 2023
    // Using the numeric API - no chrono dependency required
    let result = spa::sunrise_sunset_utc_for_horizon(
        2023,
        6,
        21,
        37.7749,   // San Francisco latitude
        -122.4194, // San Francisco longitude
        69.0,      // Delta T (seconds)
        Horizon::SunriseSunset,
    )?;

    match result {
        solar_positioning::SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        } => {
            println!("San Francisco, June 21, 2023 (UTC):");
            println!("  Sunrise:  {} hours", sunrise.hours());
            println!("  Transit:  {} hours", transit.hours());
            println!("  Sunset:   {} hours", sunset.hours());
            println!();

            // Show how to convert to day offset and hours
            let (day_offset, hours) = sunrise.day_and_hours();
            println!("Sunrise breakdown:");
            println!("  Day offset: {}", day_offset);
            println!("  Hours in day: {:.2}", hours);
        }
        solar_positioning::SunriseResult::AllDay { transit } => {
            println!("Polar day - sun never sets");
            println!("  Transit: {} hours", transit.hours());
        }
        solar_positioning::SunriseResult::AllNight { transit } => {
            println!("Polar night - sun never rises");
            println!("  Transit: {} hours", transit.hours());
        }
    }

    // Example with custom elevation angle
    println!("\nWith custom elevation angle (-1.0°):");
    let custom_result = spa::sunrise_sunset_utc(
        2023, 6, 21, 37.7749, -122.4194, 69.0, -1.0, // Custom elevation angle
    )?;

    if let solar_positioning::SunriseResult::RegularDay {
        sunrise,
        transit,
        sunset,
    } = custom_result
    {
        println!("  Sunrise:  {} hours", sunrise.hours());
        println!("  Transit:  {} hours", transit.hours());
        println!("  Sunset:   {} hours", sunset.hours());
    }

    Ok(())
}
