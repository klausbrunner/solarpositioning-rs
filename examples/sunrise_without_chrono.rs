//! Bounded event searches without chrono or heap allocation.
use solar_positioning::{Horizon, Location, SolarEvents, time::JulianDate};

fn main() -> solar_positioning::Result<()> {
    let calculator = SolarEvents::new();
    let start = JulianDate::from_utc(2023, 6, 21, 0, 0, 0.0, 69.184)?.julian_date();
    let end = start + 1.0;
    let horizon = Horizon::SunriseSunset; // Or Horizon::Custom(-4.5).
    let location = Location {
        latitude: 37.7749,
        longitude: -122.4194,
    };

    let rise = calculator.next_rise_from_julian(start, end, location, 69.184, horizon)?;
    let transit = calculator.next_transit_from_julian(start, end, location.longitude, 69.184)?;
    let set = calculator.next_set_from_julian(start, end, location, 69.184, horizon)?;
    for (name, event) in [("Rise", rise), ("Transit", transit), ("Set", set)] {
        match event {
            Some(time) => println!(
                "{name}: {:.6} hours after midnight UTC",
                (time - start) * 24.0
            ),
            None => println!("No {name} in this interval"),
        }
    }
    // To continue, use a returned event as the next search's exclusive start.
    Ok(())
}
