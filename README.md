# solar-positioning

[![CI](https://github.com/klausbrunner/solarpositioning-rs/workflows/CI/badge.svg)](https://github.com/klausbrunner/solarpositioning-rs/actions/workflows/ci.yml) [![Crates.io](https://img.shields.io/crates/v/solar-positioning?color=dodgerblue)](https://crates.io/crates/solar-positioning) [![docs.rs](https://img.shields.io/docsrs/solar-positioning)](https://docs.rs/solar-positioning)

A Rust library for finding topocentric solar coordinates, i.e. the sun's position on the sky for a given date, latitude, and longitude (and other parameters), as well as times of sunrise, sunset and twilight. Position calculations follow well-known, peer-reviewed algorithms: [SPA](https://doi.org/10.1016/j.solener.2003.12.003) by Reda and Andreas and, alternatively, [Grena/ENEA](https://doi.org/10.1016/j.solener.2012.01.024) by Grena. More than 1000 test points are included to validate against the reference code and other sources. Solar events are found by searching positions from the chosen algorithm, with SPA as the default.

> [!NOTE]
> This library is **not** based on or derived from any code published by NREL, ENEA or other parties. It implements the position algorithms as described in the respective papers.

## Status

This Rust version was originally bootstrapped from the mature [Java](https://github.com/klausbrunner/solarpositioning) project. The algorithmic core is fairly well tested; APIs may still evolve. Breaking changes may occur in minor version updates. You'll probably want to pin to a specific version in production code.

## Usage

```sh
cargo add solar-positioning
cargo add chrono # For the DateTime examples below
```

### Requirements

Rust 1.85+ (edition 2024). Minimal dependencies. Supports `std` (default) and `no_std` with `libm`.

**Feature flags:**

- `std` (default): Standard library, native math, and `alloc`
- `chrono` (default): `DateTime` API (disable for the numeric API)
- `alloc` (enabled by `std`): Collect event lists, also available with `no_std`
- `libm`: `no_std` support

### Code

`SolarPositions` and `SolarEvents` are immutable, reusable calculators, both using SPA by default.
They share `Location` for geographic coordinates.

```rust
use chrono::{DateTime, FixedOffset};
use solar_positioning::{Location, SolarPositions};

let time = "2025-06-21T12:00:00+02:00".parse::<DateTime<FixedOffset>>()?;
let location = Location { latitude: 48.21, longitude: 16.37 };
let positions = SolarPositions::new();
let position = positions.at(&time, location, 190.0, 69.0, None)?;

println!("Azimuth: {:.1}°, elevation: {:.1}°",
    position.azimuth(), position.elevation_angle());
```

Height is in metres above sea level; pass `0.0` for sea level. Delta T is in seconds.
Pass `None` for an unrefracted position, or `Some(RefractionCorrection::standard())`
for refraction under standard atmospheric conditions. Use
`RefractionCorrection::new(pressure, temperature)?` for local pressure in hPa and temperature in °C.
`SolarPositions::grena3()` selects Grena3, which requires zero height.

For many locations at one time, prepare the time-dependent calculations once:

```rust
let prepared = positions.for_time(&time, 69.0)?;
for location in locations {
    let position = prepared.at(location, 0.0, None)?;
}
```

The returned `PreparedPositions` owns its cached values and can outlive the calculator
and timestamp. Both models support this API. Calculators and prepared values are
`Copy`, `Send` and `Sync`, with no heap allocation.

Calendar constructors use the proleptic Gregorian calendar, including before 1582.
Without chrono, both models accept a `JulianDate`, which includes delta T:

```rust
use solar_positioning::{time::JulianDate, Location, SolarPositions};

let time = JulianDate::from_utc(2025, 6, 21, 12, 0, 0.0, 69.0)?;
let location = Location { latitude: 48.21, longitude: 16.37 };
let positions = SolarPositions::new();
let position = positions.at_from_julian(time, location, 190.0, None)?;
let prepared = positions.for_time_from_julian(time);
```

`SolarEvents` finds all events in a local calendar date:

```rust
use chrono::{FixedOffset, NaiveDate};
use solar_positioning::{SolarEvents, Horizon, Location, delta_t};

let date = NaiveDate::from_ymd_opt(2026, 6, 21).unwrap();
let zone = FixedOffset::east_opt(2 * 3600).unwrap();
let location = Location { latitude: 48.21, longitude: 16.37 };
let calculator = SolarEvents::new();
let day = calculator.for_date(
    date, &zone, location,
    delta_t::estimate_from_date_like(date)?, Horizon::SunriseSunset,
)?;
println!("Sunrises: {:?}", day.rises);
println!("Transits: {:?}", day.transits);
println!("Sunsets: {:?}", day.sets);
println!("Continuous daylight: {}", day.always_above());
```

Each list can be empty or contain several events. Times use the requested zone and
lie within the date, including its start and excluding the following date. Pass an
IANA zone from `chrono-tz` for daylight-saving rules, repeated dates and skipped dates.
`state_at_start` describes the Sun relative to the horizon, including `OnHorizon`
within the numerical rounding allowance. `always_above()` and `always_below()` are
false for an empty date.

Use `Horizon::CivilTwilight`, `NauticalTwilight`, `AstronomicalTwilight`, or
`Custom(elevation_in_degrees)` for other crossings. `for_date_multiple` accepts an
iterator of horizons and returns a vector of `(Horizon, Events)` pairs in input order,
sharing the transit calculation. Repeated horizons are preserved.

For an arbitrary interval, `next_rise`, `next_set` and `next_transit` return
`Result<Option<DateTime<Tz>>>` and borrow the start and end timestamps.
`None` means the event is absent. Searches exclude `start` and include `end`;
resume from a returned time to find the following event.
Rise and set are independent of transit.

```rust
let sunrise = calculator.next_rise(
    &start, &end, location, 69.184, Horizon::SunriseSunset,
)?;
```

`SolarEvents::grena3()` selects Grena3. `SolarEvents::with_provider(provider, years)`
accepts a function or closure taking a `JulianDate` and `Location` and returning
`Result<EventPosition>`, plus an inclusive supported year range. Providers supply
unrefracted sea-level elevation and local hour angle, and must satisfy the smoothness
contract documented on that constructor.

Without chrono, use `next_rise_from_julian`, `next_set_from_julian` and
`next_transit_from_julian`, which accept and return continuous UT Julian dates.
These searches need no heap allocation. With `alloc`, `for_interval_from_julian`
collects all events in a supplied half-open interval. The calendar collection
helpers require both `chrono` and `alloc` (both are enabled by default).

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

- SPA: Maximum accuracy, reference algorithm, works for historic dates
- Grena3: Simple, very fast, often accurate enough (2010-2110 CE timeframe)

Both are fast in absolute terms. The ~10× speed difference only matters for bulk calculations.

### Solar event accuracy

`SolarEvents` searches the chosen model's unrefracted, topocentric solar-centre
positions at sea level.
The standard sunrise/sunset horizon is 50 arcminutes (about 0.833°) below the
geometric horizon, allowing for average refraction and the Sun's apparent radius.
Twilight and custom horizons use their selected angle without adding refraction.

Crossing brackets and date assignment have one-millisecond resolution. This is
numerical precision, not observed-event accuracy: shallow crossings amplify the
position model's angular uncertainty, and weather, observer elevation and terrain
can shift observed sunrise by minutes ([USNO](https://aa.usno.navy.mil/faq/RST_defs)).
A tangency is not a crossing, and events less than one millisecond apart need not
be distinguished. Transit is upper meridian passage, not necessarily maximum elevation.

Event searches use continuous time. Both UT and TT must stay within the model's
supported years (-2000 through 6000 for SPA, 2010 through 2110 for Grena3).
UTC approximates UT1, and delta T is held constant within each search.

The search combines interval subdivision, [interpolation error bounds](https://dlmf.nist.gov/3.3.E5)
and [ITP refinement](https://doi.org/10.1145/3423597). Its conservative curvature estimate
accounts for daily rotation and slower solar motion. It is checked against reference data,
not formally guaranteed for every input.
[Astronomy Engine](https://github.com/cosinekitty/astronomy/blob/master/source/js/astronomy.ts)
uses related adaptive-search ideas with a speed bound instead of curvature.

### Delta T

Delta T (ΔT) is the difference between terrestrial time and UT1 ([Wikipedia](<https://en.wikipedia.org/wiki/ΔT_(timekeeping)>)). For many applications it's negligible (~70 seconds in 2025). For maximum accuracy, use observed values (available from US Naval Observatory) or estimates.

`delta_t::estimate(decimal_year)` returns an estimate in seconds. Use
`delta_t::estimate_from_date(year, month)` for a calendar month, or
`delta_t::estimate_from_date_like(date)` with chrono dates.

The estimates use polynomials originally published by [Espenak and Meeus](http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html) and [updated by Espenak in 2014](https://www.eclipsewise.com/help/deltatpoly2014.html), with custom replacement branches from 2015 onwards.
The [derivation and comparisons](https://klaus.brunners.name/posts/delta-t-polynomials/) describe the fit and its
limitations. Future values remain uncertain, and extrapolation beyond 2100 is particularly speculative.

## License

Licensed under the MIT License.
