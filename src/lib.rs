//! # Solar Positioning Library
//!
//! High-accuracy solar positioning algorithms for calculating sun position and sunrise/sunset times.

#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(docsrs, feature(doc_cfg))]
//!
//! This library provides implementations of two complementary solar positioning algorithms:
//! - **SPA** (Solar Position Algorithm): NREL's authoritative algorithm (±0.0003°, years -2000 to 6000)
//! - **Grena3**: Simplified algorithm (±0.01°, years 2010-2110, ~10x faster)
//!
//! In addition, [`delta_t`] estimates Delta T (ΔT) based on the work of F. Espenak & J. Meeus,
//! with updated fits from 2015 onwards.
//!
//! ## Features
//!
//! - Multiple configurations: `std` or `no_std`, with or without `chrono`, math via native or `libm`
//! - Maximum accuracy: Authentic NREL SPA implementation, validated against reference data
//! - Prepared calculations for many locations at one time, for both models
//! - Thread-safe: Stateless, immutable data structures
//!
//! ## Feature Flags
//!
//! - `std` (default): Native math functions and `alloc`
//! - `alloc` (enabled by `std`): Collect event lists, also available with `no_std`
//! - `chrono` (default): Enable `DateTime<Tz>` based convenience API
//! - `libm`: Use pure Rust math for `no_std` environments
//!
//! **Configuration examples:**
//! ```toml
//! # Default: std + chrono (most convenient)
//! solar-positioning = "0.7"
//!
//! # Minimal std (no chrono, smallest dependency tree)
//! solar-positioning = { version = "0.7", default-features = false, features = ["std"] }
//!
//! # no_std + chrono (embedded with DateTime support)
//! solar-positioning = { version = "0.7", default-features = false, features = ["libm", "chrono"] }
//!
//! # Minimal no_std (pure numeric API)
//! solar-positioning = { version = "0.7", default-features = false, features = ["libm"] }
//! ```
//!
//! ## References
//!
//! - Reda, I.; Andreas, A. (2003). Solar position algorithm for solar radiation applications.
//!   Solar Energy, 76(5), 577-589. DOI: <https://doi.org/10.1016/j.solener.2003.12.003>
//! - Grena, R. (2012). Five new algorithms for the computation of sun position from 2010 to 2110.
//!   Solar Energy, 86(5), 1323-1337. DOI: <https://doi.org/10.1016/j.solener.2012.01.024>
//!
//! ## Quick Start
//!
//! ### Solar positions
//! ```rust
//! # #[cfg(feature = "chrono")] {
//! use chrono::{DateTime, Utc};
//! use solar_positioning::{Location, RefractionCorrection, SolarPositions};
//!
//! let time = "2026-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
//! let location = Location { latitude: 48.21, longitude: 16.37 };
//! let positions = SolarPositions::new();
//! let position = positions.at(
//!     &time, location, 190.0, 69.0, Some(RefractionCorrection::standard()),
//! )?;
//! println!("Elevation: {:.3}°", position.elevation_angle());
//!
//! // Reuse the time-dependent calculations across locations.
//! let prepared = positions.for_time(&time, 69.0)?;
//! let unrefracted = prepared.at(location, 0.0, None)?;
//! # }
//! # Ok::<(), solar_positioning::Error>(())
//! ```
//!
//! `SolarPositions::grena3()` selects Grena3. Height is in metres (zero is required
//! for Grena3), delta T is in seconds, and `None` omits refraction. The prepared
//! value owns its cached data and needs no allocation.
//!
//! ### Positions without chrono
//! ```rust
//! use solar_positioning::{time::JulianDate, Location, SolarPositions};
//!
//! let time = JulianDate::from_utc(2026, 6, 21, 12, 0, 0.0, 69.0)?;
//! let location = Location { latitude: 48.21, longitude: 16.37 };
//! let positions = SolarPositions::new();
//! let position = positions.at_from_julian(time, location, 190.0, None)?;
//! let prepared = positions.for_time_from_julian(time);
//! assert_eq!(position, prepared.at(location, 190.0, None)?);
//! # Ok::<(), solar_positioning::Error>(())
//! ```
//!
//! ### Solar events
//! ```rust
//! # #[cfg(all(feature = "chrono", feature = "alloc"))] {
//! use chrono::{NaiveDate, Utc};
//! use solar_positioning::{SolarEvents, Horizon, Location};
//!
//! let location = Location { latitude: 78.216667, longitude: 15.633333 };
//! let day = SolarEvents::new().for_date(
//!     NaiveDate::from_ymd_opt(2020, 4, 16).unwrap(), &Utc,
//!     location, 69.184, Horizon::SunriseSunset,
//! ).unwrap();
//! assert_eq!(day.rises.len(), 2); // A date need not contain just one solar cycle.
//! # }
//! ```
//!
//! Without chrono or allocation, use continuous UT Julian dates:
//! ```rust
//! use solar_positioning::{SolarEvents, Horizon, Location};
//!
//! let location = Location { latitude: 48.21, longitude: 16.37 };
//! let rise = SolarEvents::new().next_rise_from_julian(
//!     2_460_389.5, 2_460_390.5, location, 69.184, Horizon::SunriseSunset,
//! ).unwrap();
//! assert!(rise.is_some());
//! ```
//!
//! ## Algorithms
//!
//! ### SPA (Solar Position Algorithm)
//!
//! Based on the NREL algorithm by Reda & Andreas (2003). Provides the highest accuracy
//! with uncertainties of ±0.0003 degrees, suitable for applications requiring precise
//! solar positioning over long time periods.
//!
//! ### Grena3
//!
//! A simplified algorithm optimized for years 2010-2110. Approximately 10 times faster
//! than SPA while maintaining good accuracy (maximum error 0.01°).
//!
//! ## Coordinate System
//!
//! - **Azimuth**: 0° = North, measured clockwise (0° to 360°)
//! - **Zenith angle**: 0° = directly overhead (zenith), 90° = horizon (0° to 180°)
//! - **Elevation angle**: 0° = horizon, 90° = directly overhead (-90° to +90°)

#![deny(missing_docs)]
#![deny(unsafe_code)]
#![warn(clippy::pedantic, clippy::nursery, clippy::cargo, clippy::all)]
#![allow(
    clippy::module_name_repetitions,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cargo_common_metadata,
    clippy::multiple_crate_versions, // Acceptable for dev-dependencies
    clippy::float_cmp, // Exact comparisons of mathematical constants in tests
)]

// Public API exports - core types only
pub use crate::error::{Error, Result};
#[cfg(feature = "alloc")]
pub use crate::events::Events;
pub use crate::events::{EventPosition, HorizonState, SolarEvents};
pub use crate::positions::{PreparedPositions, SolarPositions};
pub use crate::types::{Horizon, Location, RefractionCorrection, SolarPosition};

#[cfg(feature = "alloc")]
extern crate alloc;

// Algorithm modules
pub mod events;
mod grena3;
mod spa;

// Supporting modules
pub mod delta_t;
mod error;
pub mod time;
mod types;

// Internal modules
mod math;
mod positions;
