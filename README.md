# solar-positioning

[![CI](https://github.com/klausbrunner/solarpositioning-rs/workflows/CI/badge.svg)](https://github.com/klausbrunner/solarpositioning-rs/actions/workflows/ci.yml) [![Crates.io](https://img.shields.io/crates/v/solar-positioning?color=dodgerblue)](https://crates.io/crates/solar-positioning) [![docs.rs](https://img.shields.io/docsrs/solar-positioning)](https://docs.rs/solar-positioning)

A Rust library for finding topocentric solar coordinates, i.e. the sun's position on the sky for a given date, latitude, and longitude (and other parameters), as well as times of sunrise, sunset and twilight. Calculations strictly follow well-known, peer-reviewed algorithms: [SPA](http://dx.doi.org/10.1016/j.solener.2003.12.003) by Reda and Andreas and, alternatively, [Grena/ENEA](http://dx.doi.org/10.1016/j.solener.2012.01.024) by Grena. More than 1000 test points are included to validate against the reference code and other sources.

## Status

This library is a port of the Java [solarpositioning](https://github.com/klausbrunner/solarpositioning) library to idiomatic Rust and produces identical results. Unlike the very mature Java project though, it should be considered in *beta* status for now.

## Usage

```sh
cargo add solar-positioning
```

### Requirements

Rust 1.85 or newer. Minimal dependencies (chrono for date/time, thiserror for errors). Supports `no_std` with `libm` feature.

### Code

The API is intentionally "flat", comprising a handful of functions and simple structs as results.

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

For better performance when calculating positions for many coordinates at the same time, you can use the split functions to avoid recalculating time-dependent values:

```rust
use solar_positioning::spa;

let datetime = "2025-06-21T12:00:00+02:00".parse().unwrap();

// Calculate time-dependent parts once
let time_dependent = spa::spa_time_dependent_parts(datetime, 69.0).unwrap();

// Reuse for multiple coordinates
for (lat, lon) in [(48.21, 16.37), (52.52, 13.40), (45.46, 9.19)] {
    let position = spa::spa_with_time_dependent_parts(
        lat, lon, 0.0, None, &time_dependent
    ).unwrap();
    println!("{:.1}, {:.1}: {:.1}°", lat, lon, position.azimuth());
}
```

The `spa` module includes functions to calculate the times of sunrise, sun transit, and sunset in one fell swoop. The actual return type depends on the type of day (regular day, polar day, polar night).

```rust
use chrono::{DateTime, FixedOffset};
use solar_positioning::{spa, types::SunriseResult, Horizon, time::DeltaT};

// Northern location (Tromsø, Norway) in local time
let datetime = "2025-06-21T00:00:00+02:00".parse::<DateTime<FixedOffset>>().unwrap();

let result = spa::sunrise_sunset_for_horizon(
    datetime,
    69.65, // latitude
    18.96, // longitude
    DeltaT::estimate_from_date_like(datetime).unwrap(), // delta T
    Horizon::SunriseSunset
).unwrap();

match result {
    SunriseResult::RegularDay { sunrise, transit, sunset } => {
        println!("Sunrise: {}, Transit: {}, Sunset: {}", sunrise, transit, sunset);
    }
    _ => {
        println!("no sunrise or sunset today!");
    }
}
```

Twilight start and end times can be obtained like sunrise and sunset, but assuming a different horizon:

```rust
use chrono::{DateTime, FixedOffset};
use solar_positioning::{spa, Horizon, time::DeltaT};

let datetime = "2025-06-21T00:00:00+02:00".parse::<DateTime<FixedOffset>>().unwrap();

let result = spa::sunrise_sunset_for_horizon(
    datetime,
    70.978, // latitude
    25.974, // longitude
    DeltaT::estimate_from_date_like(datetime).unwrap(), // delta T
    Horizon::CivilTwilight
).unwrap();
```

See the documentation for more functions.

### Examples

The library includes several examples demonstrating different use cases:

```bash
# Basic solar position calculation
cargo run --example basic_usage

# Sunrise/sunset with different twilight types
cargo run --example sunrise_sunset

# Compare SPA vs Grena3 algorithms
cargo run --example grena3_comparison
```

### Which position algorithm should I use?

* For many applications, `grena3` should work just fine. It's simple, fast, and pretty accurate for a time window from 2010 to 2110 CE.
* If you're looking for maximum accuracy or need to calculate for historic dates, use `spa`. It's widely considered a reference algorithm for solar positioning, being very accurate and usable in a very large time window.

While `grena3` is about an order of magnitude faster than `spa`, both algorithms are very fast in absolute terms. The performance difference will only matter for bulk calculations (e.g. for many locations or times).

### Notes on sunrise, sunset, and twilight

* Calculation is based on the usual correction of 0.833° on the zenith angle, i.e. sunrise and sunset are assumed to occur when the center of the solar disc is 50 arc-minutes below the horizon. While commonly used, this fixed value fails to account for the varying effects of atmospheric refraction. Calculated and apparent sunrise and sunset times may easily differ by several minutes (cf. [Wilson 2018](https://doi.org/10.37099/mtu.dc.etdr/697)).
* As a general note on accuracy, Jean Meeus advises that "giving rising or setting times .. more accurately than to the nearest minute makes no sense" (_Astronomical Algorithms_). Errors increase the farther the position from the equator, i.e. values for polar regions are much less reliable.
* The SPA sunset/sunrise algorithm is one of the most accurate ones around. Results of this implementation correspond very closely to the [NOAA calculator](http://www.esrl.noaa.gov/gmd/grad/solcalc/)'s.

### What's this "delta T" thing?

See [Wikipedia](https://en.wikipedia.org/wiki/ΔT_(timekeeping)) for an explanation. For many simple applications, and particularly for sunrise and sunset, this value could be negligible as it's just over a minute (about 70 seconds as of 2025). However, if you're looking for maximum accuracy, you should use an observed value (available from e.g. the US Naval Observatory) or at least a solid estimate.

The `time::DeltaT` type provides an estimator based on polynomials fitting a number of observed (or extrapolated) historical values, published by [Espenak and Meeus](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html) in 2007 and slightly updated by [Espenak](https://www.eclipsewise.com/help/deltatpoly2014.html) in 2014.

As of September 2025, extrapolated values from this estimator are slightly high (about 2 seconds). This gap will widen in the coming decades (cf. [Morrison et al. 2021](https://royalsocietypublishing.org/doi/10.1098/rspa.2020.0776)). Still, the estimates should work sufficiently well for most applications.

### Is it thread-safe?

Yes. All functions are stateless and safe to call concurrently.

## License

Licensed under the MIT License.
