//! Core data types for solar positioning calculations.

use crate::Result;
use crate::error::{check_azimuth, check_pressure, check_temperature, check_zenith_angle};

const SUNRISE_SUNSET_ELEVATION: f64 = -50.0 / 60.0;

/// Geographic coordinates in degrees, with longitude positive east.
///
/// Event calculations use sea level. Calculators validate latitude in [-90, 90]
/// and longitude in [-180, 180], rejecting non-finite coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Location {
    /// Latitude in degrees, positive north.
    pub latitude: f64,
    /// Longitude in degrees, positive east.
    pub longitude: f64,
}

/// Predefined elevation angles for sunrise/sunset calculations.
///
/// Corresponds to different twilight definitions for consistent sunrise, sunset, and twilight calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Horizon {
    /// Standard sunrise/sunset (sun's upper limb touches horizon, accounting for refraction)
    SunriseSunset,
    /// Civil twilight (sun is 6° below horizon)
    CivilTwilight,
    /// Nautical twilight (sun is 12° below horizon)
    NauticalTwilight,
    /// Astronomical twilight (sun is 18° below horizon)
    AstronomicalTwilight,
    /// Custom elevation angle
    Custom(f64),
}

impl Horizon {
    /// Gets the elevation angle in degrees for this horizon definition.
    ///
    /// Negative values indicate the sun is below the horizon.
    #[must_use]
    pub const fn elevation_angle(&self) -> f64 {
        match self {
            Self::SunriseSunset => SUNRISE_SUNSET_ELEVATION, // Accounts for refraction and sun's radius
            Self::CivilTwilight => -6.0,
            Self::NauticalTwilight => -12.0,
            Self::AstronomicalTwilight => -18.0,
            Self::Custom(angle) => *angle,
        }
    }
}

/// Atmospheric conditions for refraction correction in solar position calculations.
///
/// Atmospheric refraction bends light rays, causing the apparent sun position to differ
/// from its true geometric position by up to ~0.6° near the horizon.
///
/// # Example
/// ```
/// # use solar_positioning::RefractionCorrection;
/// // Standard atmospheric conditions at sea level
/// let standard = RefractionCorrection::standard();
/// assert_eq!(standard.pressure(), 1013.25);
/// assert_eq!(standard.temperature(), 15.0);
///
/// // Custom conditions for high altitude or different climate
/// let custom = RefractionCorrection::new(900.0, -5.0).unwrap();
/// assert_eq!(custom.pressure(), 900.0);
/// assert_eq!(custom.temperature(), -5.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefractionCorrection {
    /// Atmospheric pressure in millibars (hPa)
    pressure: f64,
    /// Temperature in degrees Celsius
    temperature: f64,
}

impl RefractionCorrection {
    /// Creates a new refraction correction with the specified atmospheric conditions.
    ///
    /// # Errors
    /// Returns `InvalidPressure` unless pressure is finite, positive and at most 2000 hPa.
    /// Returns `InvalidTemperature` unless temperature is finite, greater than -273°C and at most 100°C.
    ///
    /// # Example
    /// ```
    /// # use solar_positioning::RefractionCorrection;
    /// let correction = RefractionCorrection::new(1013.25, 15.0).unwrap();
    /// assert_eq!(correction.pressure(), 1013.25);
    /// assert_eq!(correction.temperature(), 15.0);
    /// ```
    pub fn new(pressure: f64, temperature: f64) -> Result<Self> {
        check_pressure(pressure)?;
        check_temperature(temperature)?;
        Ok(Self {
            pressure,
            temperature,
        })
    }

    /// Creates refraction correction using standard atmospheric conditions.
    ///
    /// Uses standard sea-level conditions:
    /// - Pressure: 1013.25 millibars (standard atmosphere)
    /// - Temperature: 15.0°C (59°F)
    ///
    /// # Example
    /// ```
    /// # use solar_positioning::RefractionCorrection;
    /// let standard = RefractionCorrection::standard();
    /// assert_eq!(standard.pressure(), 1013.25);
    /// assert_eq!(standard.temperature(), 15.0);
    /// ```
    #[must_use]
    pub const fn standard() -> Self {
        Self {
            pressure: 1013.25,
            temperature: 15.0,
        }
    }

    /// Gets the atmospheric pressure in millibars.
    #[must_use]
    pub const fn pressure(&self) -> f64 {
        self.pressure
    }

    /// Gets the temperature in degrees Celsius.
    #[must_use]
    pub const fn temperature(&self) -> f64 {
        self.temperature
    }
}

/// Solar position in topocentric coordinates.
///
/// Represents the sun's position as seen from a specific point on Earth's surface.
/// Uses the standard astronomical coordinate system where:
/// - Azimuth: 0° = North, measured clockwise to 360°
/// - Zenith angle: 0° = directly overhead (zenith), 90° = horizon, 180° = nadir
/// - Elevation angle: 90° = directly overhead, 0° = horizon, -90° = nadir
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarPosition {
    /// Azimuth angle in degrees (0° to 360°, 0° = North, increasing clockwise)
    azimuth: f64,
    /// Zenith angle in degrees (0° to 180°, 0° = zenith, 90° = horizon)
    zenith_angle: f64,
}

impl SolarPosition {
    /// Creates a new solar position from azimuth and zenith angle.
    ///
    /// # Errors
    /// Returns error if azimuth or zenith angles are outside valid ranges.
    ///
    /// # Example
    /// ```
    /// # use solar_positioning::SolarPosition;
    /// let position = SolarPosition::new(180.0, 30.0).unwrap();
    /// assert_eq!(position.azimuth(), 180.0);
    /// assert_eq!(position.zenith_angle(), 30.0);
    /// assert_eq!(position.elevation_angle(), 60.0);
    /// ```
    pub fn new(azimuth: f64, zenith_angle: f64) -> Result<Self> {
        let normalized_azimuth = check_azimuth(azimuth)?;
        let validated_zenith = check_zenith_angle(zenith_angle)?;

        Ok(Self {
            azimuth: normalized_azimuth,
            zenith_angle: validated_zenith,
        })
    }

    /// Gets the azimuth angle in degrees (0° to 360°, 0° = North, increasing clockwise).
    #[must_use]
    pub const fn azimuth(&self) -> f64 {
        self.azimuth
    }

    /// Gets the zenith angle in degrees (0° to 180°, 0° = zenith, 90° = horizon).
    #[must_use]
    pub const fn zenith_angle(&self) -> f64 {
        self.zenith_angle
    }

    /// Gets the elevation angle in degrees.
    ///
    /// This is the complement of the zenith angle: elevation = 90° - zenith.
    #[must_use]
    pub fn elevation_angle(&self) -> f64 {
        90.0 - self.zenith_angle
    }

    /// Checks if the sun is above the horizon (elevation angle > 0°).
    #[must_use]
    pub fn is_sun_up(&self) -> bool {
        self.elevation_angle() > 0.0
    }

    /// Checks if the sun is at or below the horizon (elevation angle ≤ 0°).
    #[must_use]
    pub fn is_sun_down(&self) -> bool {
        self.elevation_angle() <= 0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_horizon_elevation_angles() {
        assert_eq!(Horizon::SunriseSunset.elevation_angle(), -50.0 / 60.0);
        assert_eq!(Horizon::CivilTwilight.elevation_angle(), -6.0);
        assert_eq!(Horizon::NauticalTwilight.elevation_angle(), -12.0);
        assert_eq!(Horizon::AstronomicalTwilight.elevation_angle(), -18.0);

        assert_eq!(Horizon::Custom(-3.0).elevation_angle(), -3.0);
    }

    #[test]
    fn test_solar_position_creation() {
        let pos = SolarPosition::new(180.0, 45.0).unwrap();
        assert_eq!(pos.azimuth(), 180.0);
        assert_eq!(pos.zenith_angle(), 45.0);
        assert_eq!(pos.elevation_angle(), 45.0);
        assert!(pos.is_sun_up());
        assert!(!pos.is_sun_down());

        // Test normalization
        let pos = SolarPosition::new(-90.0, 90.0).unwrap();
        assert_eq!(pos.azimuth(), 270.0);
        assert_eq!(pos.elevation_angle(), 0.0);

        // Test validation
        assert!(SolarPosition::new(0.0, -1.0).is_err());
        assert!(SolarPosition::new(0.0, 181.0).is_err());
    }

    #[test]
    fn test_solar_position_sun_state() {
        let above_horizon = SolarPosition::new(180.0, 30.0).unwrap();
        assert!(above_horizon.is_sun_up());
        assert!(!above_horizon.is_sun_down());

        let on_horizon = SolarPosition::new(180.0, 90.0).unwrap();
        assert!(!on_horizon.is_sun_up());
        assert!(on_horizon.is_sun_down());

        let below_horizon = SolarPosition::new(180.0, 120.0).unwrap();
        assert!(!below_horizon.is_sun_up());
        assert!(below_horizon.is_sun_down());
    }

    #[test]
    fn test_refraction_correction() {
        // Test standard conditions
        let standard = RefractionCorrection::standard();
        assert_eq!(standard.pressure(), 1013.25);
        assert_eq!(standard.temperature(), 15.0);

        // Test custom conditions
        let custom = RefractionCorrection::new(1000.0, 20.0).unwrap();
        assert_eq!(custom.pressure(), 1000.0);
        assert_eq!(custom.temperature(), 20.0);

        // Test validation
        assert!(RefractionCorrection::new(-1.0, 15.0).is_err()); // Invalid pressure
        assert!(RefractionCorrection::new(1013.25, -300.0).is_err()); // Invalid temperature
        assert!(RefractionCorrection::new(3000.0, 15.0).is_err()); // Too high pressure
        assert!(RefractionCorrection::new(1013.25, 150.0).is_err()); // Too high temperature
    }
}
