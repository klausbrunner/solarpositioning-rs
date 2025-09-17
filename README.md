# solar-positioning

[![CI](https://github.com/klausbrunner/solarpositioning-rs/workflows/CI/badge.svg)](https://github.com/klausbrunner/solarpositioning-rs/actions/workflows/ci.yml) [![Crates.io](https://img.shields.io/crates/v/solar-positioning?color=dodgerblue)](https://crates.io/crates/solar-positioning) [![docs.rs](https://img.shields.io/docsrs/solar-positioning)](https://docs.rs/solar-positioning)

A Rust library for finding topocentric solar coordinates, i.e. the sun's position on the sky for a given date, latitude, and longitude (and other parameters), as well as times of sunrise and sunset. Calculations strictly follow well-known, peer-reviewed algorithms: [SPA](http://dx.doi.org/10.1016/j.solener.2003.12.003) by Reda and Andreas and, alternatively, [Grena/ENEA](http://dx.doi.org/10.1016/j.solener.2012.01.024) by Grena. More than 1000 test points are included to validate against the reference code and other sources.

## Status

This library is a port of the Java [solarpositioning](https://github.com/klausbrunner/solarpositioning) library to idiomatic Rust and should produce identical results. Unlike the very mature Java project though, it should be considered in *beta* status for now.

## Usage

### Cargo.toml

```toml
[dependencies]
solar-positioning = "0.2.1"
```

### Requirements

Rust 1.85 or newer. Minimal dependencies (chrono for date/time, thiserror for errors).

### Code

The API is intentionally "flat", comprising a handful of functions and simple structs as results.
To get refraction-corrected topocentric coordinates:

```rust
use chrono::{DateTime, FixedOffset};
use solar_positioning::{spa, time::DeltaT};

// Use timezone-aware datetime (Vienna, Central European Time)
let datetime = "2025-06-21T12:00:00+02:00".parse::<DateTime<FixedOffset>>().unwrap();

// replace spa with grena3 as needed
let position = spa::solar_position(
    datetime,
    48.21,   // latitude (degrees)
    16.37,   // longitude (degrees)
    190.0,   // elevation (m)
    DeltaT::estimate_from_date_like(&datetime).unwrap(), // delta T (s, ~70s for 2025)
    1010.0,  // avg. air pressure (hPa)
    11.0     // avg. air temperature (°C)
).unwrap();

println!("{:?}", position);
```

The spa module includes functions to calculate the times of sunrise, sun transit, and sunset in one fell swoop. The actual return type depends on the type of day (regular day, polar day, polar night).

```rust
use chrono::{DateTime, FixedOffset};
use solar_positioning::{spa, types::SunriseResult, Horizon, time::DeltaT};

// Arctic location (Tromsø, Norway) in local time
let datetime = "2025-06-21T00:00:00+02:00".parse::<DateTime<FixedOffset>>().unwrap();

let result = spa::sunrise_sunset_for_horizon(
    datetime,
    70.978, // latitude
    25.974, // longitude
    DeltaT::estimate_from_date_like(&datetime).unwrap(), // delta T
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
    DeltaT::estimate_from_date_like(&datetime).unwrap(), // delta T
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

* For many applications, Grena3 should work just fine. It's simple, fast, and pretty accurate for a time window from 2010 to 2110 CE.
* If you're looking for maximum accuracy or need to calculate for historic dates, use SPA. It's widely considered a reference algorithm for solar positioning, being very accurate and usable in a very large time window.

While Grena3 is about an order of magnitude faster than SPA, both algorithms are very fast in absolute terms.

### Notes on sunrise, sunset, and twilight

* Calculation is based on the usual correction of 0.833° on the zenith angle, i.e. sunrise and sunset are assumed to occur when the center of the solar disc is 50 arc-minutes below the horizon. While commonly used, this fixed value fails to account for the varying effects of atmospheric refraction. Calculated and apparent sunrise and sunset times may easily differ by several minutes (cf. [Wilson 2018](https://doi.org/10.37099/mtu.dc.etdr/697)).
* As a general note on accuracy, Jean Meeus advises that "giving rising or setting times .. more accurately than to the nearest minute makes no sense" (_Astronomical Algorithms_). Errors increase the farther the position from the equator, i.e. values for polar regions are much less reliable.
* The SPA sunset/sunrise algorithm is one of the most accurate ones around. Results of this implementation correspond very closely to the [NOAA calculator](http://www.esrl.noaa.gov/gmd/grad/solcalc/)'s.

### What's this "delta T" thing?

See [Wikipedia](https://en.wikipedia.org/wiki/ΔT_(timekeeping)) for an explanation. For many simple applications, and particularly for sunrise and sunset, this value could be negligible as it's just over a minute (about 70 seconds as of 2025). However, if you're looking for maximum accuracy, you should use an observed value (available from e.g. the US Naval Observatory) or at least a solid estimate.

The DeltaT module provides an estimator based on polynomials fitting a number of observed (or extrapolated) historical values, published by [Espenak and Meeus](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html) in 2007 and slightly updated by [Espenak](https://www.eclipsewise.com/help/deltatpoly2014.html) in 2014.

As of September 2025, extrapolated values from this estimator are slightly high (about 2 seconds). This gap will widen in the coming decades (cf. [Morrison et al. 2021](https://royalsocietypublishing.org/doi/10.1098/rspa.2020.0776)). Still, the estimates should work sufficiently well for most applications.

### Is it thread-safe?

Yes. All functions are stateless and safe to call concurrently.

## License

Licensed under the MIT License.
