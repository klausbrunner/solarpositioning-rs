# solar-positioning

[![CI](https://github.com/klausbrunner/solarpositioning-rs/workflows/CI/badge.svg)](https://github.com/klausbrunner/solarpositioning-rs/actions/workflows/ci.yml) [![Crates.io](https://img.shields.io/crates/v/solar-positioning?color=dodgerblue)](https://crates.io/crates/solar-positioning) [![docs.rs](https://img.shields.io/docsrs/solar-positioning)](https://docs.rs/solar-positioning)

A Rust library for finding topocentric solar coordinates, i.e. the sun's position on the sky for a given date, latitude, and longitude (and other parameters), as well as times of sunrise, sunset and twilight. Calculations strictly follow well-known, peer-reviewed algorithms: [SPA](http://dx.doi.org/10.1016/j.solener.2003.12.003) by Reda and Andreas and, alternatively, [Grena/ENEA](http://dx.doi.org/10.1016/j.solener.2012.01.024) by Grena. More than 1000 test points are included to validate against the reference code and other sources.

> [!NOTE]
> This library is **not** based on or derived from code published by NREL, ENEA or other parties. It is an implementation precisely following the algorithms described in the respective papers.

## Status

While the core algorithms are stable and well-tested, the API is still evolving. Breaking changes may occur in minor version updates. Please pin to a specific version in production code.

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
    let pos = spa::spa_with_time_dependent_parts(&time_dependent, lat, lon, 0.0, None).unwrap();
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

For twilight, use `Horizon::CivilTwilight`, `Horizon::NauticalTwilight`, or `Horizon::AstronomicalTwilight`.

### Examples

```bash
cargo run --example basic_usage          # Solar position
cargo run --example sunrise_sunset       # Sunrise/sunset/twilight
cargo run --example grena3_comparison    # SPA vs Grena3
```

### Which algorithm?

- `spa`: Maximum accuracy, reference algorithm, works for historic dates
- `grena3`: Simple, very fast, often accurate enough (2010-2110 CE timeframe)

Both are fast in absolute terms. The ~10× speed difference only matters for bulk calculations.

### Sunrise/sunset accuracy notes

- Uses standard 0.833° correction (solar disc 50 arc-minutes below horizon). Atmospheric refraction varies, so calculated times may differ from observed by several minutes ([Wilson 2018](https://doi.org/10.37099/mtu.dc.etdr/697)).
- Jean Meeus advises giving times "more accurately than to the nearest minute makes no sense". Errors increase toward poles.
- Results match the [NOAA calculator](http://www.esrl.noaa.gov/gmd/grad/solcalc/) closely.

#### Divergence from the NREL SPA reference code

The library follows the procedure in the SPA paper: sidereal time is evaluated at 0 **UT** (A.2.1)
while the geocentric α/δ for sunrise/sunset interpolation are evaluated at 0 **TT** for D−1/D/D+1 (A.2.2). The NREL
reference code (`spa.c`) resets ΔT to zero when building those intermediate ephemerides, effectively keeping
them in UT. This Rust code preserves the supplied ΔT to stay faithful to the published algorithm rather than
the C code. As a consequence, sunrise/sunset times differ slightly from `spa.c` but should line up better with
high-precision ephemerides (JPL Horizons, USNO almanacs, etc.).

### Delta T

Delta T (ΔT) is the difference between terrestrial time and UT1 ([Wikipedia](https://en.wikipedia.org/wiki/ΔT_(timekeeping))). For many applications it's negligible (~70 seconds in 2025). For maximum accuracy, use observed values (available from US Naval Observatory) or estimates.

The `time::DeltaT` estimator uses polynomial fits from [Espenak and Meeus](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html) (2007, updated 2014). Current extrapolated values are slightly high (~2 seconds). This gap will widen ([Morrison et al. 2021](https://royalsocietypublishing.org/doi/10.1098/rspa.2020.0776)). However, this should not matter for most applications.

## License

Licensed under the MIT License.
