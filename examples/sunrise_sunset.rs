//! Independent rise/set events and twilight for a local calendar date.
use jiff::{civil::date, tz::TimeZone};
use solar_positioning::{Horizon, Location, SolarEvents};

fn main() -> solar_positioning::Result<()> {
    let calculator = SolarEvents::new(); // Or SolarEvents::grena3().
    let date = date(2020, 4, 16);
    let zone = TimeZone::get("Arctic/Longyearbyen").unwrap();
    let location = Location {
        latitude: 78.216667,
        longitude: 15.633333,
    };
    let horizons = [
        Horizon::SunriseSunset,
        Horizon::CivilTwilight,
        Horizon::NauticalTwilight,
        Horizon::AstronomicalTwilight,
    ];

    for (horizon, day) in calculator.for_date_multiple(date, &zone, location, 69.184, horizons)? {
        println!("{date}, {horizon:?}");
        for time in day.rises {
            println!("  Rise:    {time}");
        }
        for time in day.transits {
            println!("  Transit: {time}");
        }
        for time in day.sets {
            println!("  Set:     {time}");
        }
    }
    Ok(())
}
