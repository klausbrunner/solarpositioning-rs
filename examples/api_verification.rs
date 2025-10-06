//! Verify the public API surface makes sense.
//!
//! This example demonstrates all the main ways users can calculate sunrise/sunset.

use solar_positioning::{spa, Horizon, HoursUtc, SunriseResult};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== API Verification ===\n");

    // 1. Non-chrono API with explicit elevation angle
    println!("1. Non-chrono API with explicit elevation angle:");
    let result1: SunriseResult<HoursUtc> = spa::sunrise_sunset_utc(
        2023, 6, 21, 37.7749, -122.4194, 69.0, -0.833, // explicit angle
    )?;
    print_hours_result(&result1);

    // 2. Non-chrono API with Horizon enum (convenience)
    println!("\n2. Non-chrono API with Horizon enum:");
    let result2: SunriseResult<HoursUtc> = spa::sunrise_sunset_utc_for_horizon(
        2023,
        6,
        21,
        37.7749,
        -122.4194,
        69.0,
        Horizon::SunriseSunset,
    )?;
    print_hours_result(&result2);

    // 3. Different horizons
    println!("\n3. Different horizon types:");
    for horizon in [
        Horizon::SunriseSunset,
        Horizon::CivilTwilight,
        Horizon::NauticalTwilight,
        Horizon::AstronomicalTwilight,
    ] {
        let result =
            spa::sunrise_sunset_utc_for_horizon(2023, 6, 21, 37.7749, -122.4194, 69.0, horizon)?;
        if let SunriseResult::RegularDay {
            sunrise,
            transit: _,
            sunset,
        } = result
        {
            println!(
                "  {:?}: {:.2} - {:.2} hours",
                horizon,
                sunrise.hours(),
                sunset.hours()
            );
        }
    }

    // 4. Working with HoursUtc
    println!("\n4. Working with HoursUtc:");
    if let SunriseResult::RegularDay {
        sunrise,
        transit,
        sunset,
    } = result1
    {
        println!("  Sunrise raw hours: {}", sunrise.hours());
        let (day, hours) = sunrise.day_and_hours();
        println!("  Sunrise normalized: day={}, hours={:.2}", day, hours);

        println!("  Transit raw hours: {}", transit.hours());
        println!("  Sunset raw hours: {}", sunset.hours());
    }

    // 5. Chrono API (if feature enabled)
    #[cfg(feature = "chrono")]
    {
        use chrono::{FixedOffset, NaiveDate, TimeZone};

        println!("\n5. Chrono API:");
        let date = FixedOffset::east_opt(-7 * 3600)
            .unwrap()
            .from_local_datetime(
                &NaiveDate::from_ymd_opt(2023, 6, 21)
                    .unwrap()
                    .and_hms_opt(0, 0, 0)
                    .unwrap(),
            )
            .unwrap();

        let result_chrono = spa::sunrise_sunset_for_horizon(
            date,
            37.7749,
            -122.4194,
            69.0,
            Horizon::SunriseSunset,
        )?;

        if let SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        } = result_chrono
        {
            println!("  Sunrise: {}", sunrise);
            println!("  Transit: {}", transit);
            println!("  Sunset:  {}", sunset);
        }
    }

    println!("\n✓ API verification complete!");
    Ok(())
}

fn print_hours_result(result: &SunriseResult<HoursUtc>) {
    match result {
        SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        } => {
            println!("  Sunrise: {:.2} hours", sunrise.hours());
            println!("  Transit: {:.2} hours", transit.hours());
            println!("  Sunset:  {:.2} hours", sunset.hours());
        }
        SunriseResult::AllDay { transit } => {
            println!("  Polar day (transit: {:.2} hours)", transit.hours());
        }
        SunriseResult::AllNight { transit } => {
            println!("  Polar night (transit: {:.2} hours)", transit.hours());
        }
    }
}
