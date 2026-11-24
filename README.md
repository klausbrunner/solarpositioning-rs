# solar-positioning

[![CI](https://github.com/klausbrunner/solarpositioning-rs/workflows/CI/badge.svg)](https://github.com/klausbrunner/solarpositioning-rs/actions/workflows/ci.yml) [![Crates.io](https://img.shields.io/crates/v/solar-positioning?color=dodgerblue)](https://crates.io/crates/solar-positioning) [![docs.rs](https://img.shields.io/docsrs/solar-positioning)](https://docs.rs/solar-positioning)

A Rust library for topocentric solar positions and sunrise, sunset, transit and
twilight times. It implements [SPA](https://doi.org/10.1016/j.solener.2003.12.003)
by Reda and Andreas and [Grena3](https://doi.org/10.1016/j.solener.2012.01.024) by Grena.
Tests compare SPA positions against NREL reference data and solar events against
JPL and USNO data.

Event times use a [separate numerical search](https://klaus.brunners.name/posts/spa-sunrise-revisited/)
over positions from the selected model, not SPA's built-in event calculation.
It combines interval subdivision with ITP root refinement.

> [!NOTE]
> This library is **not** based on or derived from any code published by NREL, ENEA or other parties. It implements the position algorithms as described in the respective papers.

## Features

- SPA (default): ±0.0003° accuracy, years -2000 through 6000.
- Grena3: ±0.01° accuracy, years 2010 through 2110, roughly 10× faster.
- Timezone-aware and numeric Julian-date APIs, including `no_std` support.
- Optional observer height and atmospheric refraction, and elevation-crossing searches.
- Prepared calculations for many locations at one time. Positions and next-event
  searches need no heap allocation.

## Usage

Install the library and Jiff for the timestamp and calendar examples:

```sh
cargo add solar-positioning jiff@1
```

Requires Rust 1.85+. See the [API reference](https://docs.rs/solar-positioning)
for method contracts and additional configurations.

- `std` (default): Standard library, native math, and `alloc`
- `jiff` (default): Jiff timestamps and, with `alloc`, local calendar events
- `alloc` (enabled by `std`): Collect event lists, also available with `no_std`
- `libm`: `no_std` support

### Solar position

```rust
use jiff::ZonedDateTime;
use solar_positioning::{Location, SolarPositions, delta_t};

let zoned = "2025-06-21T12:00:00+02:00[Europe/Vienna]".parse::<ZonedDateTime>()?;
let time = zoned.timestamp();
let location = Location { latitude: 48.21, longitude: 16.37 };
let positions = SolarPositions::new();
let deltat = delta_t::estimate_from_civil_date(zoned.date())?;
let position = positions.at(time, location, 190.0, deltat, None)?;

println!("Azimuth: {:.1}°, elevation: {:.1}°",
    position.azimuth(), position.elevation_angle());
```

Height is in metres above sea level, delta T (TT minus UT1) in seconds, and angles in degrees.
`None` omits refraction; use `Some(RefractionCorrection::standard())` to include it.
`SolarPositions::grena3()` selects Grena3.

### Many locations at one time

Reuse the time-dependent work from the position example:

```rust
let locations = [location, Location { latitude: 52.52, longitude: 13.41 }];
let prepared = positions.for_time(time, deltat)?;
for location in locations {
    let position = prepared.at(location, 0.0, None)?;
}
```

### Sunrise, sunset and twilight

```rust
use jiff::{civil::date, tz::TimeZone};
use solar_positioning::{SolarEvents, Horizon, Location, delta_t};

let date = date(2026, 6, 21);
let zone = TimeZone::get("Europe/Vienna")?;
let location = Location { latitude: 48.21, longitude: 16.37 };
let calculator = SolarEvents::new();
let day = calculator.for_date(
    date, &zone, location,
    delta_t::estimate_from_civil_date(date)?, Horizon::SunriseSunset,
)?;
println!("Sunrises: {:?}", day.rises);
println!("Transits: {:?}", day.transits);
println!("Sunsets: {:?}", day.sets);
println!("Continuous daylight: {}", day.always_above());
```

Lists may be empty or contain multiple events. Use a twilight `Horizon` or
`Horizon::Custom(angle)` for other crossings, and `for_date_multiple` for several
horizons. Calendar events are returned in the supplied time zone. Jiff handles
the time zone database; this library does not bundle one.
`SolarEvents::grena3()` selects Grena3.

### Next elevation crossing

Using the time and location from the position example, find the next rising
crossing of an apparent solar-centre elevation:

```rust
use jiff::SignedDuration;
use solar_positioning::{ElevationTarget, RefractionCorrection};

let end = time + SignedDuration::from_hours(24);
let target = ElevationTarget {
    angle: 10.0,
    height: 190.0,
    refraction: Some(RefractionCorrection::standard()),
};
let rise = positions.next_elevation_rise(time, end, location, deltat, target)?;
```

Angles are apparent with refraction, geometric without it. The result is `None` if
the interval contains no crossing; `next_elevation_set` finds the descending crossing.
Use `.to_zoned_datetime(zone)` on a returned timestamp for local time.

### Without Jiff

Disable default features and enable `std`, or `libm` for `no_std`:

```rust
use solar_positioning::{time::JulianDate, Location, SolarPositions};

let time = JulianDate::from_utc(2025, 6, 21, 12, 0, 0.0, 69.0)?;
let location = Location { latitude: 48.21, longitude: 16.37 };
let position = SolarPositions::new().at_from_julian(time, location, 190.0, None)?;
```

Event searches also have `*_from_julian` methods.

## Limitations

Supported years apply to both UT and TT. Calendar inputs use the proleptic
Gregorian calendar; UTC approximates UT1. Jiff parses a leap second (`:60`) as
`:59`; numeric UTC components reject `:60`. Grena3 requires zero observer height.
Height changes parallax, not the visible horizon.

Standard events use sea-level solar-centre positions, with sunrise and sunset including
the conventional refraction and solar-radius allowance. Search resolution is 100 ms.
Weather and terrain can shift observed sunrise by minutes ([USNO](https://aa.usno.navy.mil/faq/RST_defs)).
Transit is meridian passage, not necessarily maximum elevation.

Delta T estimates follow [Espenak and Meeus](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html)
with [updated fits](https://klaus.brunners.name/posts/delta-t-polynomials/).
Future estimates are (obviously) uncertain, particularly beyond 2100.

## Examples

```bash
cargo run --example basic_usage          # Solar position
cargo run --example sunrise_sunset       # Sunrise/sunset/twilight
cargo run --example grena3_comparison    # SPA vs Grena3
cargo run --example sunrise_from_julian --no-default-features --features std
cargo run --example no_std_usage --no-default-features --features libm
```

## Status

Originally bootstrapped from the mature [Java](https://github.com/klausbrunner/solarpositioning)
project. APIs may still change in minor releases; pin a version for production use.

## License

Licensed under the MIT License.
