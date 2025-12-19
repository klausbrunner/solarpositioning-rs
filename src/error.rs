//! Error types for the solar positioning library.

use crate::math::normalize_degrees_0_to_360;
use core::fmt;

/// Result type alias for operations in this crate.
pub type Result<T> = core::result::Result<T, Error>;

/// Errors that can occur during solar position calculations.
#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    /// Invalid latitude value (must be between -90 and +90 degrees).
    InvalidLatitude {
        /// The invalid latitude value provided.
        value: f64,
    },
    /// Invalid longitude value (must be between -180 and +180 degrees).
    InvalidLongitude {
        /// The invalid longitude value provided.
        value: f64,
    },
    /// Invalid elevation angle for sunrise/sunset calculations.
    InvalidElevationAngle {
        /// The invalid elevation angle value provided.
        value: f64,
    },
    /// Invalid pressure value for atmospheric refraction calculations.
    InvalidPressure {
        /// The invalid pressure value provided.
        value: f64,
    },
    /// Invalid temperature value for atmospheric refraction calculations.
    InvalidTemperature {
        /// The invalid temperature value provided.
        value: f64,
    },
    /// Invalid date/time for the algorithm's valid range.
    InvalidDateTime {
        /// Description of the date/time constraint violation.
        message: &'static str,
    },
    /// Numerical computation error (e.g., convergence failure).
    ComputationError {
        /// Description of the computation error.
        message: &'static str,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLatitude { value } => {
                write!(
                    f,
                    "invalid latitude {value}° (must be between -90° and +90°)"
                )
            }
            Self::InvalidLongitude { value } => {
                write!(
                    f,
                    "invalid longitude {value}° (must be between -180° and +180°)"
                )
            }
            Self::InvalidElevationAngle { value } => {
                write!(
                    f,
                    "invalid elevation angle {value}° (must be between -90° and +90°)"
                )
            }
            Self::InvalidPressure { value } => {
                write!(f, "invalid pressure {value} mbar (must be positive)")
            }
            Self::InvalidTemperature { value } => {
                write!(
                    f,
                    "invalid temperature {value}°C (must be above absolute zero)"
                )
            }
            Self::InvalidDateTime { message } => {
                write!(f, "invalid date/time: {message}")
            }
            Self::ComputationError { message } => {
                write!(f, "computation error: {message}")
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {}

impl Error {
    /// Creates an invalid latitude error.
    #[must_use]
    pub const fn invalid_latitude(value: f64) -> Self {
        Self::InvalidLatitude { value }
    }

    /// Creates an invalid longitude error.
    #[must_use]
    pub const fn invalid_longitude(value: f64) -> Self {
        Self::InvalidLongitude { value }
    }

    /// Creates an invalid elevation angle error.
    #[must_use]
    pub const fn invalid_elevation_angle(value: f64) -> Self {
        Self::InvalidElevationAngle { value }
    }

    /// Creates an invalid pressure error.
    #[must_use]
    pub const fn invalid_pressure(value: f64) -> Self {
        Self::InvalidPressure { value }
    }

    /// Creates an invalid temperature error.
    #[must_use]
    pub const fn invalid_temperature(value: f64) -> Self {
        Self::InvalidTemperature { value }
    }

    /// Creates an invalid date/time error.
    #[must_use]
    pub const fn invalid_datetime(message: &'static str) -> Self {
        Self::InvalidDateTime { message }
    }

    /// Creates a computation error.
    #[must_use]
    pub const fn computation_error(message: &'static str) -> Self {
        Self::ComputationError { message }
    }
}

/// Validates latitude is within the valid range (-90 to +90 degrees).
///
/// # Errors
/// Returns `InvalidLatitude` if latitude is outside -90 to +90 degrees.
pub fn check_latitude(latitude: f64) -> Result<()> {
    if !(-90.0..=90.0).contains(&latitude) {
        return Err(Error::invalid_latitude(latitude));
    }
    Ok(())
}

/// Validates longitude is within the valid range (-180 to +180 degrees).
///
/// # Errors
/// Returns `InvalidLongitude` if longitude is outside -180 to +180 degrees.
pub fn check_longitude(longitude: f64) -> Result<()> {
    if !(-180.0..=180.0).contains(&longitude) {
        return Err(Error::invalid_longitude(longitude));
    }
    Ok(())
}

/// Validates both latitude and longitude are within valid ranges.
///
/// # Errors
/// Returns `InvalidLatitude` or `InvalidLongitude` for out-of-range coordinates.
pub fn check_coordinates(latitude: f64, longitude: f64) -> Result<()> {
    check_latitude(latitude)?;
    check_longitude(longitude)?;
    Ok(())
}

/// Validates pressure is positive and reasonable for atmospheric calculations.
///
/// # Errors
/// Returns `InvalidPressure` if pressure is not between 1 and 2000 hPa.
pub fn check_pressure(pressure: f64) -> Result<()> {
    if !pressure.is_finite() || pressure <= 0.0 || pressure > 2000.0 {
        return Err(Error::invalid_pressure(pressure));
    }
    Ok(())
}

/// Validates temperature is above absolute zero and reasonable for atmospheric calculations.
///
/// # Errors
/// Returns `InvalidTemperature` if temperature is outside -273.15 to 100°C.
pub fn check_temperature(temperature: f64) -> Result<()> {
    if !(-273.15..=100.0).contains(&temperature) {
        return Err(Error::invalid_temperature(temperature));
    }
    Ok(())
}

/// Validates and normalizes an azimuth angle to the range [0, 360) degrees.
///
/// # Errors
/// Returns `ComputationError` if azimuth is not finite.
pub fn check_azimuth(azimuth: f64) -> Result<f64> {
    if !azimuth.is_finite() {
        return Err(Error::computation_error("azimuth is not finite"));
    }
    Ok(normalize_degrees_0_to_360(azimuth))
}

/// Validates a zenith angle to be within the range [0, 180] degrees.
///
/// # Errors
/// Returns `ComputationError` if zenith angle is not finite or outside valid range.
pub fn check_zenith_angle(zenith: f64) -> Result<f64> {
    if !zenith.is_finite() {
        return Err(Error::computation_error("zenith angle is not finite"));
    }
    if !(0.0..=180.0).contains(&zenith) {
        return Err(Error::computation_error(
            "zenith angle must be between 0° and 180°",
        ));
    }
    Ok(zenith)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_latitude_validation() {
        assert!(check_latitude(0.0).is_ok());
        assert!(check_latitude(90.0).is_ok());
        assert!(check_latitude(-90.0).is_ok());
        assert!(check_latitude(45.5).is_ok());

        assert!(check_latitude(91.0).is_err());
        assert!(check_latitude(-91.0).is_err());
        assert!(check_latitude(f64::NAN).is_err());
        assert!(check_latitude(f64::INFINITY).is_err());
    }

    #[test]
    fn test_longitude_validation() {
        assert!(check_longitude(0.0).is_ok());
        assert!(check_longitude(180.0).is_ok());
        assert!(check_longitude(-180.0).is_ok());
        assert!(check_longitude(122.5).is_ok());

        assert!(check_longitude(181.0).is_err());
        assert!(check_longitude(-181.0).is_err());
        assert!(check_longitude(f64::NAN).is_err());
        assert!(check_longitude(f64::INFINITY).is_err());
    }

    #[test]
    fn test_pressure_validation() {
        assert!(check_pressure(1013.25).is_ok());
        assert!(check_pressure(1000.0).is_ok());
        assert!(check_pressure(500.0).is_ok());

        assert!(check_pressure(0.0).is_err());
        assert!(check_pressure(-100.0).is_err());
        assert!(check_pressure(3000.0).is_err());
        assert!(check_pressure(f64::NAN).is_err());
        assert!(check_pressure(f64::INFINITY).is_err());
    }

    #[test]
    fn test_temperature_validation() {
        assert!(check_temperature(15.0).is_ok());
        assert!(check_temperature(0.0).is_ok());
        assert!(check_temperature(-40.0).is_ok());
        assert!(check_temperature(50.0).is_ok());

        assert!(check_temperature(-300.0).is_err());
        assert!(check_temperature(150.0).is_err());
    }

    #[test]
    #[cfg(feature = "std")]
    fn test_error_display() {
        let err = Error::invalid_latitude(95.0);
        assert_eq!(
            err.to_string(),
            "invalid latitude 95° (must be between -90° and +90°)"
        );

        let err = Error::invalid_longitude(185.0);
        assert_eq!(
            err.to_string(),
            "invalid longitude 185° (must be between -180° and +180°)"
        );

        let err = Error::computation_error("convergence failed");
        assert_eq!(err.to_string(), "computation error: convergence failed");
    }

    #[test]
    fn test_check_azimuth() {
        assert!(check_azimuth(0.0).is_ok());
        assert!(check_azimuth(180.0).is_ok());
        assert!(check_azimuth(360.0).is_ok());
        assert!(check_azimuth(450.0).is_ok());
        assert!(check_azimuth(-90.0).is_ok());

        // Check normalization
        assert_eq!(check_azimuth(-90.0).unwrap(), 270.0);
        assert_eq!(check_azimuth(450.0).unwrap(), 90.0);

        assert!(check_azimuth(f64::NAN).is_err());
        assert!(check_azimuth(f64::INFINITY).is_err());
    }

    #[test]
    fn test_check_zenith_angle() {
        assert!(check_zenith_angle(0.0).is_ok());
        assert!(check_zenith_angle(90.0).is_ok());
        assert!(check_zenith_angle(180.0).is_ok());

        assert!(check_zenith_angle(-1.0).is_err());
        assert!(check_zenith_angle(181.0).is_err());
        assert!(check_zenith_angle(f64::NAN).is_err());
        assert!(check_zenith_angle(f64::INFINITY).is_err());
    }
}
