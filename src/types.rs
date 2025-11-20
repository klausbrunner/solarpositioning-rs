//! Core data types for solar positioning calculations.

use crate::error::{check_azimuth, check_pressure, check_temperature, check_zenith_angle};
use crate::math::floor;
use crate::{Error, Result};

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
            Self::SunriseSunset => -0.83337, // Accounts for refraction and sun's radius
            Self::CivilTwilight => -6.0,
            Self::NauticalTwilight => -12.0,
            Self::AstronomicalTwilight => -18.0,
            Self::Custom(angle) => *angle,
        }
    }

    /// Creates a custom horizon with the specified elevation angle.
    ///
    /// # Errors
    /// Returns `InvalidElevationAngle` if elevation is outside -90 to +90 degrees.
    pub fn custom(elevation_degrees: f64) -> Result<Self> {
        if !(-90.0..=90.0).contains(&elevation_degrees) {
            return Err(Error::invalid_elevation_angle(elevation_degrees));
        }
        Ok(Self::Custom(elevation_degrees))
    }
}

impl Eq for Horizon {}

impl core::hash::Hash for Horizon {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        match self {
            Self::SunriseSunset => 0.hash(state),
            Self::CivilTwilight => 1.hash(state),
            Self::NauticalTwilight => 2.hash(state),
            Self::AstronomicalTwilight => 3.hash(state),
            Self::Custom(angle) => {
                4.hash(state);
                // Normalize -0.0 and +0.0 so hashing remains consistent with PartialEq
                let normalized = if *angle == 0.0 { 0.0 } else { *angle };
                normalized.to_bits().hash(state);
            }
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
/// # use solar_positioning::types::RefractionCorrection;
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
    /// Returns `InvalidPressure` or `InvalidTemperature` for out-of-range values.
    ///
    /// # Example
    /// ```
    /// # use solar_positioning::types::RefractionCorrection;
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
    /// # use solar_positioning::types::RefractionCorrection;
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
    /// # use solar_positioning::types::SolarPosition;
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

/// Hours since midnight UTC that can extend beyond a single day.
///
/// Used for sunrise/sunset times without the chrono dependency.
/// Values represent hours since midnight UTC (0 UT) for the calculation date:
/// - Negative values indicate the previous day
/// - 0.0 to < 24.0 indicates the current day
/// - ≥ 24.0 indicates the next day
///
/// # Example
/// ```
/// # use solar_positioning::types::HoursUtc;
/// let morning = HoursUtc::from_hours(6.5); // 06:30 current day
/// let late_evening = HoursUtc::from_hours(23.5); // 23:30 current day
/// let after_midnight = HoursUtc::from_hours(24.5); // 00:30 next day
/// let before_midnight_prev = HoursUtc::from_hours(-0.5); // 23:30 previous day
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HoursUtc(f64);

impl HoursUtc {
    /// Creates a new `HoursUtc` from hours since midnight UTC.
    ///
    /// Values can be negative (previous day) or ≥ 24.0 (next day).
    #[must_use]
    pub const fn from_hours(hours: f64) -> Self {
        Self(hours)
    }

    /// Gets the raw hours value.
    ///
    /// Can be negative (previous day) or ≥ 24.0 (next day).
    #[must_use]
    pub const fn hours(&self) -> f64 {
        self.0
    }

    /// Gets the day offset and normalized hours (0.0 to < 24.0).
    ///
    /// # Returns
    /// Tuple of (`day_offset`, `hours_in_day`) where:
    /// - `day_offset`: whole days offset from the calculation date (negative = previous days, positive = following days)
    /// - `hours_in_day`: 0.0 to < 24.0
    ///
    /// # Example
    /// ```
    /// # use solar_positioning::types::HoursUtc;
    /// let time = HoursUtc::from_hours(25.5);
    /// let (day_offset, hours) = time.day_and_hours();
    /// assert_eq!(day_offset, 1);
    /// assert!((hours - 1.5).abs() < 1e-10);
    /// ```
    #[must_use]
    pub fn day_and_hours(&self) -> (i32, f64) {
        let hours = self.0;
        if !hours.is_finite() {
            return (0, hours);
        }

        let mut day_offset_raw = floor(hours / 24.0);
        let mut normalized_hours = hours - day_offset_raw * 24.0;

        if normalized_hours < 0.0 {
            normalized_hours += 24.0;
            day_offset_raw -= 1.0;
        } else if normalized_hours >= 24.0 {
            normalized_hours -= 24.0;
            day_offset_raw += 1.0;
        }

        let day_offset = if day_offset_raw >= f64::from(i32::MAX) {
            i32::MAX
        } else if day_offset_raw <= f64::from(i32::MIN) {
            i32::MIN
        } else {
            day_offset_raw as i32
        };

        (day_offset, normalized_hours)
    }
}

/// Result of sunrise/sunset calculations for a given day.
///
/// Solar events can vary significantly based on location and time of year,
/// especially at extreme latitudes where polar days and nights occur.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(
    feature = "std",
    doc = "Default generic parameter is `()`; chrono helpers return `SunriseResult<chrono::DateTime<Tz>>`."
)]
pub enum SunriseResult<T = ()> {
    /// Regular day with distinct sunrise, transit (noon), and sunset times
    RegularDay {
        /// Time of sunrise
        sunrise: T,
        /// Time of solar transit (when sun crosses meridian, solar noon)
        transit: T,
        /// Time of sunset
        sunset: T,
    },
    /// Polar day - sun remains above the specified horizon all day
    AllDay {
        /// Time of solar transit (closest approach to zenith)
        transit: T,
    },
    /// Polar night - sun remains below the specified horizon all day
    AllNight {
        /// Time of solar transit (when sun is highest, though still below horizon)
        transit: T,
    },
}

impl<T> SunriseResult<T> {
    /// Gets the transit time (solar noon) for any sunrise result.
    pub const fn transit(&self) -> &T {
        match self {
            Self::RegularDay { transit, .. }
            | Self::AllDay { transit }
            | Self::AllNight { transit } => transit,
        }
    }

    /// Checks if this represents a regular day with sunrise and sunset.
    pub const fn is_regular_day(&self) -> bool {
        matches!(self, Self::RegularDay { .. })
    }

    /// Checks if this represents a polar day (sun never sets).
    pub const fn is_polar_day(&self) -> bool {
        matches!(self, Self::AllDay { .. })
    }

    /// Checks if this represents a polar night (sun never rises).
    pub const fn is_polar_night(&self) -> bool {
        matches!(self, Self::AllNight { .. })
    }

    /// Gets sunrise time if this is a regular day.
    pub const fn sunrise(&self) -> Option<&T> {
        if let Self::RegularDay { sunrise, .. } = self {
            Some(sunrise)
        } else {
            None
        }
    }

    /// Gets sunset time if this is a regular day.
    pub const fn sunset(&self) -> Option<&T> {
        if let Self::RegularDay { sunset, .. } = self {
            Some(sunset)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_horizon_elevation_angles() {
        assert_eq!(Horizon::SunriseSunset.elevation_angle(), -0.83337);
        assert_eq!(Horizon::CivilTwilight.elevation_angle(), -6.0);
        assert_eq!(Horizon::NauticalTwilight.elevation_angle(), -12.0);
        assert_eq!(Horizon::AstronomicalTwilight.elevation_angle(), -18.0);

        let custom = Horizon::custom(-3.0).unwrap();
        assert_eq!(custom.elevation_angle(), -3.0);

        assert!(Horizon::custom(-95.0).is_err());
        assert!(Horizon::custom(95.0).is_err());
    }

    #[test]
    #[cfg(feature = "std")]
    fn test_horizon_hash_normalizes_zero_sign() {
        use std::collections::HashSet;

        let mut set = HashSet::new();
        set.insert(Horizon::Custom(0.0));
        set.insert(Horizon::Custom(-0.0));

        assert_eq!(set.len(), 1, "hashing should treat +0.0 and -0.0 equally");
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
    fn test_sunrise_result_regular_day() {
        use chrono::{DateTime, Utc};

        let sunrise = "2023-06-21T05:30:00Z".parse::<DateTime<Utc>>().unwrap();
        let transit = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
        let sunset = "2023-06-21T18:30:00Z".parse::<DateTime<Utc>>().unwrap();

        let result = SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        };

        assert!(result.is_regular_day());
        assert!(!result.is_polar_day());
        assert!(!result.is_polar_night());
        assert_eq!(result.transit(), &transit);
        assert_eq!(result.sunrise(), Some(&sunrise));
        assert_eq!(result.sunset(), Some(&sunset));
    }

    #[test]
    fn test_sunrise_result_polar_day() {
        use chrono::{DateTime, Utc};

        let transit = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
        let result = SunriseResult::AllDay { transit };

        assert!(!result.is_regular_day());
        assert!(result.is_polar_day());
        assert!(!result.is_polar_night());
        assert_eq!(result.transit(), &transit);
        assert_eq!(result.sunrise(), None);
        assert_eq!(result.sunset(), None);
    }

    #[test]
    fn test_sunrise_result_polar_night() {
        use chrono::{DateTime, Utc};

        let transit = "2023-12-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
        let result = SunriseResult::AllNight { transit };

        assert!(!result.is_regular_day());
        assert!(!result.is_polar_day());
        assert!(result.is_polar_night());
        assert_eq!(result.transit(), &transit);
        assert_eq!(result.sunrise(), None);
        assert_eq!(result.sunset(), None);
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
