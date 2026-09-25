# solar-positioning

[![CI](https://github.com/klausbrunner/solarpositioning-rs/workflows/CI/badge.svg)](https://github.com/klausbrunner/solarpositioning-rs/actions/workflows/ci.yml) [![Crates.io](https://img.shields.io/crates/v/solar-positioning?color=dodgerblue)](https://crates.io/crates/solar-positioning) [![docs.rs](https://img.shields.io/docsrs/solar-positioning)](https://docs.rs/solar-positioning)

A Rust library for finding topocentric solar coordinates, i.e. the sun's position on the sky for a given date, latitude, and longitude (and other parameters), as well as times of sunrise, sunset and twilight. Calculations strictly follow well-known, peer-reviewed algorithms: [SPA](http://dx.doi.org/10.1016/j.solener.2003.12.003) by Reda and Andreas and, alternatively, [Grena/ENEA](http://dx.doi.org/10.1016/j.solener.2012.01.024) by Grena. More than 1000 test points are included to validate against the reference code and other sources.

> [!NOTE]
> This library is **not** based on or derived from code published by NREL, ENEA or other parties. It implements the algorithms described in the respective papers, with minimal adjustments documented below.

## Status

This Rust version was originally bootstrapped from the mature [Java](https://github.com/klausbrunner/solarpositioning) project. The algorithmic core is fairly well tested; APIs may still evolve. Breaking changes may occur in minor version updates. You'll probably want to pin to a specific version in production code.

## Usage

```sh
cargo add solar-positioning
```

### Requirements

Rust 1.70+. Minimal dependencies. Supports `std` (default) and `no_std` with `libm`.

**Feature flags:**

- `std` (default): Standard library, native math
- `chrono` (default): `DateTime` API (disable for pure numeric `JulianDate` API)
- `libm`: `no_std` support

### Code

Functions are organized by algorithm (`spa` or `grena3` modules). Results use simple structs and enums.

```rust
use chrono::{DateTime, FixedOffset};
use solar_positioning::spa;

let datetime = "2025-06-21T12:00:00+02:00".parse::<DateTime<FixedOffset>>().unwrap();

let position = spa::solar_position(
    datetime,
    48.21,   // latitude
    16.37,   // longitude
    190.0,   // elevation (m)
    69.0,    // delta T (seconds, ~70 for 2025)
    None     // no atmospheric refraction
).unwrap();

println!("Azimuth: {:.1}°, Elevation: {:.1}°",
    position.azimuth(), position.elevation_angle());
```

Without `chrono`, use the numeric `JulianDate` API:

```rust
use solar_positioning::{spa, time::JulianDate, RefractionCorrection};

let jd = JulianDate::from_utc(2025, 6, 21, 12, 0, 0.0, 69.0).unwrap();
let position = spa::solar_position_from_julian(
    jd, 48.21, 16.37, 190.0, Some(RefractionCorrection::standard())
).unwrap();
```

For multiple coordinates at the same time, calculate time-dependent parts once (SPA only):

```rust
let time_dependent = spa::spa_time_dependent_parts(datetime, 69.0).unwrap();
for (lat, lon) in [(48.21, 16.37), (52.52, 13.40)] {
    let pos = spa::spa_with_time_dependent_parts(lat, lon, 0.0, None, &time_dependent).unwrap();
}
```

Calculate sunrise, transit, and sunset (return type depends on day type: regular/polar day/polar night):

```rust
use solar_positioning::{spa, types::SunriseResult, Horizon, time::DeltaT};

let datetime = "2025-06-21T00:00:00+02:00".parse().unwrap();
let result = spa::sunrise_sunset_for_horizon(
    datetime, 69.65, 18.96,
    DeltaT::estimate_from_date_like(datetime).unwrap(),
    Horizon::SunriseSunset
).unwrap();

match result {
    SunriseResult::RegularDay { sunrise, transit, sunset } => { /* ... */ }
    _ => { /* polar day/night */ }
}
```

Returned event timestamps are in the same timezone as the input `DateTime`, but can fall on the
previous/next local calendar date when events occur near midnight (e.g., at timezone boundaries or
for twilights).
The chrono API selects the transit closest to 12:00 on the requested date's local clock
(earlier on a tie). Transit can fall on an adjacent date; the result describes one solar cycle,
not all events in a civil day.

For twilight, use `Horizon::CivilTwilight`, `Horizon::NauticalTwilight`, or `Horizon::AstronomicalTwilight`.

### Examples

```bash
cargo run --example basic_usage          # Solar position
cargo run --example sunrise_sunset       # Sunrise/sunset/twilight
cargo run --example grena3_comparison    # SPA vs Grena3
cargo run --example sunrise_without_chrono --no-default-features --features std
cargo run --example no_std_usage --no-default-features --features libm
```

The first three examples use the default `std + chrono` feature set. The last two use the numeric API.

### Which algorithm?

- `spa`: Maximum accuracy, reference algorithm, works for historic dates
- `grena3`: Simple, very fast, often accurate enough (2010-2110 CE timeframe)

Both are fast in absolute terms. The ~10× speed difference only matters for bulk calculations.

### Sunrise/sunset accuracy notes

- Sunrise and sunset use the standard solar-centre elevation of −0.833° (50 arcminutes below the geometric horizon), accounting for average atmospheric refraction and the Sun's apparent radius.
- Atmospheric variability limits the accuracy of predicted observed sunrise/sunset times: differences of a minute or more are possible, especially where the Sun crosses the horizon at a shallow angle ([USNO](https://aa.usno.navy.mil/faq/RST_defs)).
- SPA's sunrise/sunset and twilight calculations become less reliable near seasonal transitions where the Sun barely crosses the selected horizon.
- Days with only a rising or setting event are not reliably supported by SPA.

#### Difference in SPA day wrapping

Unlike SPA Appendix A.2.7, this library retains sunrise and sunset estimates’ day offsets around the selected transit instead of wrapping them independently into [0, 1). This avoids using the wrong day’s solar coordinates. SPA’s interpolation and correction equations are unchanged.

### Delta T

Delta T (ΔT) is the difference between terrestrial time and UT1 ([Wikipedia](<https://en.wikipedia.org/wiki/ΔT_(timekeeping)>)). For many applications it's negligible (~70 seconds in 2025). For maximum accuracy, use observed values (available from US Naval Observatory) or estimates.

`time::DeltaT` provides such estimates based on polynomials originally published by [Espenak and Meeus](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html) and [updated by Espenak in 2014](https://www.eclipsewise.com/help/deltatpoly2014.html), with custom replacement branches from 2015 onwards.
The [derivation and comparisons](https://klaus.brunners.name/posts/delta-t-polynomials/) describe the fit and its
limitations. Future values remain uncertain, and extrapolation beyond 2100 is particularly speculative.

## License

Licensed under the MIT License.
