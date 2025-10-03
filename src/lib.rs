//! # Solar Positioning Library
//!
//! High-accuracy solar positioning algorithms for calculating sun position and sunrise/sunset times.

#![cfg_attr(not(feature = "std"), no_std)]
//!
//! This library provides implementations of two complementary solar positioning algorithms:
//! - **SPA** (Solar Position Algorithm): NREL's authoritative algorithm (±0.0003°, years -2000 to 6000)
//! - **Grena3**: Simplified algorithm (±0.01°, years 2010-2110, ~10x faster)
//!
//! In addition, it provides an estimator for Delta T (ΔT) values based on the work of F. Espenak & J. Meeus.
//!
//! ## Features
//!
//! - Multiple configurations: `std` or `no_std`, with or without `chrono`, math via native or `libm`
//! - Maximum accuracy: Authentic NREL SPA implementation, validated against reference data
//! - Performance optimized: Split functions for bulk calculations (SPA only)
//! - Thread-safe: Stateless, immutable data structures
//!
//! ## Feature Flags
//!
//! - `std` (default): Use standard library for native math functions (usually faster than `libm`)
//! - `chrono` (default): Enable `DateTime<Tz>` based convenience API
//! - `libm`: Use pure Rust math for `no_std` environments
//!
//! **Configuration examples:**
//! ```toml
//! # Default: std + chrono (most convenient)
//! solar-positioning = "0.3"
//!
//! # Minimal std (no chrono, smallest dependency tree)
//! solar-positioning = { version = "0.3", default-features = false, features = ["std"] }
//!
//! # no_std + chrono (embedded with DateTime support)
//! solar-positioning = { version = "0.3", default-features = false, features = ["libm", "chrono"] }
//!
//! # Minimal no_std (pure numeric API)
//! solar-positioning = { version = "0.3", default-features = false, features = ["libm"] }
//! ```
//!
//! ## References
//!
//! - Reda, I.; Andreas, A. (2003). Solar position algorithm for solar radiation applications.
//!   Solar Energy, 76(5), 577-589. DOI: <http://dx.doi.org/10.1016/j.solener.2003.12.003>
//! - Grena, R. (2012). Five new algorithms for the computation of sun position from 2010 to 2110.
//!   Solar Energy, 86(5), 1323-1337. DOI: <http://dx.doi.org/10.1016/j.solener.2012.01.024>
//!
//! ## Quick Start
//!
//! ### Solar Position (with chrono)
//! ```rust
//! # #[cfg(feature = "chrono")] {
//! use solar_positioning::{spa, RefractionCorrection, time::DeltaT};
//! use chrono::{DateTime, FixedOffset};
//!
//! // Calculate sun position for Vienna at noon
//! let datetime = "2026-06-21T12:00:00+02:00".parse::<DateTime<FixedOffset>>().unwrap();
//! let position = spa::solar_position(
//!     datetime,
//!     48.21,   // Vienna latitude
//!     16.37,   // Vienna longitude
//!     190.0,   // elevation (meters)
//!     DeltaT::estimate_from_date_like(datetime).unwrap(), // delta T
//!     Some(RefractionCorrection::standard())
//! ).unwrap();
//!
//! println!("Azimuth: {:.3}°", position.azimuth());
//! println!("Elevation: {:.3}°", position.elevation_angle());
//! # }
//! ```
//!
//! ### Solar Position (numeric API, no chrono)
//! ```rust
//! use solar_positioning::{spa, time::JulianDate, RefractionCorrection};
//!
//! // Create Julian date from UTC components (2026-06-21 12:00:00 UTC + 69s ΔT)
//! let jd = JulianDate::from_utc(2026, 6, 21, 12, 0, 0.0, 69.0).unwrap();
//!
//! // Calculate sun position (works in both std and no_std)
//! let position = spa::solar_position_from_julian(
//!     jd,
//!     48.21,   // Vienna latitude
//!     16.37,   // Vienna longitude
//!     190.0,   // elevation (meters)
//!     Some(RefractionCorrection::standard())
//! ).unwrap();
//!
//! println!("Azimuth: {:.3}°", position.azimuth());
//! println!("Elevation: {:.3}°", position.elevation_angle());
//! ```
//!
//! ### Sunrise and Sunset (requires chrono)
//! ```rust
//! # #[cfg(feature = "chrono")] {
//! use solar_positioning::{spa, Horizon, time::DeltaT};
//! use chrono::{DateTime, FixedOffset};
//!
//! // Calculate sunrise/sunset for San Francisco
//! let date = "2026-06-21T00:00:00-07:00".parse::<DateTime<FixedOffset>>().unwrap();
//! let result = spa::sunrise_sunset_for_horizon(
//!     date,
//!     37.7749,  // San Francisco latitude
//!     -122.4194, // San Francisco longitude
//!     DeltaT::estimate_from_date_like(date).unwrap(),
//!     Horizon::SunriseSunset
//! ).unwrap();
//!
//! match result {
//!     solar_positioning::SunriseResult::RegularDay { sunrise, transit, sunset } => {
//!         println!("Sunrise: {}", sunrise);
//!         println!("Solar noon: {}", transit);
//!         println!("Sunset: {}", sunset);
//!     }
//!     _ => println!("No sunrise/sunset (polar day/night)"),
//! }
//! # }
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

// Public API exports
pub use crate::error::{Error, Result};
#[cfg(feature = "chrono")]
pub use crate::spa::spa_time_dependent_parts;
pub use crate::spa::{SpaTimeDependent, spa_with_time_dependent_parts};
pub use crate::types::{Horizon, RefractionCorrection, SolarPosition, SunriseResult};

// Algorithm modules
pub mod grena3;
pub mod spa;

// Core modules
pub mod error;
pub mod types;

// Internal modules
mod math;

// Public modules
pub mod time;

#[cfg(all(test, feature = "chrono"))]
mod tests {
    use super::*;
    use chrono::{DateTime, FixedOffset, TimeZone, Utc};

    #[test]
    fn test_basic_spa_calculation() {
        // Test with different timezone types
        let datetime_fixed = "2023-06-21T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let datetime_utc = Utc.with_ymd_and_hms(2023, 6, 21, 19, 0, 0).unwrap();

        let position1 = spa::solar_position(
            datetime_fixed,
            37.7749,
            -122.4194,
            0.0,
            69.0,
            Some(RefractionCorrection::standard()),
        )
        .unwrap();
        let position2 = spa::solar_position(
            datetime_utc,
            37.7749,
            -122.4194,
            0.0,
            69.0,
            Some(RefractionCorrection::standard()),
        )
        .unwrap();

        // Both should produce identical results
        assert!((position1.azimuth() - position2.azimuth()).abs() < 1e-10);
        assert!((position1.zenith_angle() - position2.zenith_angle()).abs() < 1e-10);

        assert!(position1.azimuth() >= 0.0);
        assert!(position1.azimuth() <= 360.0);
        assert!(position1.zenith_angle() >= 0.0);
        assert!(position1.zenith_angle() <= 180.0);
    }

    #[test]
    fn test_basic_grena3_calculation() {
        use chrono::{DateTime, FixedOffset, TimeZone, Utc};

        let datetime_fixed = "2023-06-21T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let datetime_utc = Utc.with_ymd_and_hms(2023, 6, 21, 19, 0, 0).unwrap();

        let position1 = grena3::solar_position(
            datetime_fixed,
            37.7749,
            -122.4194,
            69.0,
            Some(RefractionCorrection::new(1013.25, 15.0).unwrap()),
        )
        .unwrap();

        let position2 = grena3::solar_position(
            datetime_utc,
            37.7749,
            -122.4194,
            69.0,
            Some(RefractionCorrection::new(1013.25, 15.0).unwrap()),
        )
        .unwrap();

        // Both should produce identical results
        assert!((position1.azimuth() - position2.azimuth()).abs() < 1e-6);
        assert!((position1.zenith_angle() - position2.zenith_angle()).abs() < 1e-6);

        assert!(position1.azimuth() >= 0.0);
        assert!(position1.azimuth() <= 360.0);
        assert!(position1.zenith_angle() >= 0.0);
        assert!(position1.zenith_angle() <= 180.0);
    }
}
