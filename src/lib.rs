//! # Solar Positioning Library
//!
//! High-accuracy solar positioning algorithms for calculating sun position and sunrise/sunset times.
//!
//! This library provides implementations of two complementary solar positioning algorithms:
//! - **SPA** (Solar Position Algorithm): NREL's high-accuracy algorithm (±0.0003° uncertainty, years -2000 to 6000)
//! - **Grena3**: Simplified algorithm (±0.01° accuracy, years 2010-2110, ~10x faster)

//!
//! ## Features
//!
//! - Thread-safe, immutable data structures
//! - Comprehensive test suite with reference data validation
//! - Optional serialization support with Serde
//!
//! ## Quick Start
//!
//! ```rust
//! use solar_positioning::{spa, time::JulianDate, types::SolarPosition};
//! use chrono::{DateTime, FixedOffset, Utc, TimeZone};
//!
//! // Example with time calculations
//! let jd = JulianDate::from_utc(2023, 6, 21, 12, 0, 0.0, 69.0).unwrap();
//! println!("Julian Date: {:.6}", jd.julian_date());
//! println!("Julian Century: {:.6}", jd.julian_century());
//!
//! // Example with flexible timezone support - any TimeZone trait implementor
//! let datetime_fixed = "2023-06-21T12:00:00-07:00".parse::<DateTime<FixedOffset>>().unwrap();
//! let datetime_utc = Utc.with_ymd_and_hms(2023, 6, 21, 19, 0, 0).unwrap(); // Same moment
//!
//! // Both calls produce identical results
//! let position = spa::solar_position(datetime_fixed, 37.7749, -122.4194, 0.0, 69.0, 1013.25, 15.0).unwrap();
//! println!("Azimuth: {:.3}°", position.azimuth());
//! println!("Elevation: {:.3}°", position.elevation_angle());
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
//! than SPA while maintaining good accuracy (maximum error 0.01 degrees).
//!
//! ## Coordinate System
//!
//! - **Azimuth**: 0° = North, measured clockwise (0° to 360°)
//! - **Zenith angle**: 0° = directly overhead (zenith), 90° = horizon (0° to 180°)
//! - **Elevation angle**: 0° = horizon, 90° = directly overhead (-90° to 90°)

#![cfg_attr(not(feature = "std"), no_std)]
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
pub use crate::types::{Horizon, SolarPosition, SunriseResult};

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_spa_calculation() {
        use chrono::{DateTime, FixedOffset, TimeZone, Utc};

        // Test with different timezone types
        let datetime_fixed = "2023-06-21T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let datetime_utc = Utc.with_ymd_and_hms(2023, 6, 21, 19, 0, 0).unwrap();

        let position1 =
            spa::solar_position(datetime_fixed, 37.7749, -122.4194, 0.0, 69.0, 1013.25, 15.0)
                .unwrap();
        let position2 =
            spa::solar_position(datetime_utc, 37.7749, -122.4194, 0.0, 69.0, 1013.25, 15.0)
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

        let position1 = grena3::solar_position_with_refraction(
            datetime_fixed,
            37.7749,
            -122.4194,
            69.0,
            Some(1013.25),
            Some(15.0),
        )
        .unwrap();

        let position2 = grena3::solar_position_with_refraction(
            datetime_utc,
            37.7749,
            -122.4194,
            69.0,
            Some(1013.25),
            Some(15.0),
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
