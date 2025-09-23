//! Grena3 solar position algorithm implementation.
//!
//! This follows the no. 3 algorithm described by Roberto Grena of ENEA (2012).
//!
//! The algorithm is designed for the years 2010 to 2110, with a maximum error of 0.01 degrees.
//! It's approximately 10x faster than the SPA algorithm but with reduced accuracy and time range.
//!
//! Reference: Grena, R. (2012). Five new algorithms for the computation of sun position from 2010 to 2110.
//! Solar Energy, 86(5), 1323-1337. DOI: <http://dx.doi.org/10.1016/j.solener.2012.01.024>

#![allow(clippy::unreadable_literal)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::suboptimal_flops)]

use crate::error::check_coordinates;
use crate::math::{
    PI, asin, atan2, cos, degrees_to_radians, normalize_degrees_0_to_360, radians_to_degrees, sin,
    sqrt, tan,
};
use crate::{RefractionCorrection, Result, SolarPosition};
use chrono::{DateTime, Datelike, TimeZone, Timelike};

/// Calculate solar position using the Grena3 algorithm.
///
/// This is a simplified algorithm designed for years 2010-2110 with maximum error of 0.01°.
/// It's much faster than SPA but less accurate and has a limited time range.
///
/// # Arguments
/// * `datetime` - Timezone-aware date and time
/// * `latitude` - Observer latitude in degrees (-90° to +90°)
/// * `longitude` - Observer longitude in degrees (-180° to +180°)
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `refraction` - Optional atmospheric refraction correction
///
/// # Returns
/// Returns `Ok(SolarPosition)` with azimuth and zenith angles on success.
///
/// # Errors
/// Returns error for invalid coordinates (latitude outside ±90°, longitude outside ±180°).
///
/// # Panics
/// This function does not panic.
///
/// # Example
/// ```rust
/// use solar_positioning::{grena3, RefractionCorrection};
/// use chrono::{DateTime, FixedOffset};
///
/// let datetime = "2023-06-21T12:00:00-07:00".parse::<DateTime<FixedOffset>>().unwrap();
///
/// // Without refraction correction
/// let position = grena3::solar_position(
///     datetime,
///     37.7749,     // San Francisco latitude
///     -122.4194,   // San Francisco longitude
///     69.0,        // deltaT (seconds)
///     None,        // no atmospheric refraction correction
/// ).unwrap();
///
/// // With refraction correction
/// let position_with_refraction = grena3::solar_position(
///     datetime,
///     37.7749,     // San Francisco latitude
///     -122.4194,   // San Francisco longitude
///     69.0,        // deltaT (seconds)
///     Some(RefractionCorrection::standard()), // standard atmospheric conditions
/// ).unwrap();
///
/// println!("Azimuth: {:.3}°", position.azimuth());
/// println!("Elevation: {:.3}°", position.elevation_angle());
/// ```
#[allow(clippy::needless_pass_by_value)]
pub fn solar_position<Tz: TimeZone>(
    datetime: DateTime<Tz>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    refraction: Option<RefractionCorrection>,
) -> Result<SolarPosition> {
    check_coordinates(latitude, longitude)?;

    // Calculate t (days since 2000-01-01 12:00:00 TT)
    let t = calc_t(&datetime);
    let t_e = t + 1.1574e-5 * delta_t;
    let omega_at_e = 0.0172019715 * t_e;

    // Calculate apparent sun longitude (lambda)
    let lambda = -1.388803
        + 1.720279216e-2 * t_e
        + 3.3366e-2 * sin(omega_at_e - 0.06172)
        + 3.53e-4 * sin(2.0 * omega_at_e - 0.1163);

    // Calculate obliquity of ecliptic (epsilon)
    let epsilon = 4.089567e-1 - 6.19e-9 * t_e;

    let s_lambda = sin(lambda);
    let c_lambda = cos(lambda);
    let s_epsilon = sin(epsilon);
    let c_epsilon = sqrt(1.0 - s_epsilon * s_epsilon);

    // Calculate right ascension (alpha)
    let mut alpha = atan2(s_lambda * c_epsilon, c_lambda);
    if alpha < 0.0 {
        alpha += 2.0 * PI;
    }

    // Calculate declination (delta)
    let delta = asin(s_lambda * s_epsilon);

    // Calculate hour angle (H)
    let mut h = 1.7528311 + 6.300388099 * t + degrees_to_radians(longitude) - alpha;
    h = ((h + PI) % (2.0 * PI)) - PI;
    if h < -PI {
        h += 2.0 * PI;
    }

    // Calculate topocentric coordinates
    let s_phi = sin(degrees_to_radians(latitude));
    let c_phi = sqrt(1.0 - s_phi * s_phi);
    let s_delta = sin(delta);
    let c_delta = sqrt(1.0 - s_delta * s_delta);
    let s_h = sin(h);
    let c_h = cos(h);

    let s_epsilon0 = s_phi * s_delta + c_phi * c_delta * c_h;
    let e_p = asin(s_epsilon0) - 4.26e-5 * sqrt(1.0 - s_epsilon0 * s_epsilon0);
    let gamma = atan2(s_h, c_h * s_phi - s_delta * c_phi / c_delta);

    // Apply refraction correction if provided and sun is visible
    let delta_re = refraction.map_or(0.0, |correction| {
        if e_p > 0.0 {
            let pressure = correction.pressure();
            let temperature = correction.temperature();
            (0.08422 * (pressure / 1000.0))
                / ((273.0 + temperature) * tan(e_p + 0.003138 / (e_p + 0.08919)))
        } else {
            0.0
        }
    });

    let z = PI / 2.0 - e_p - delta_re;

    let azimuth = normalize_degrees_0_to_360(radians_to_degrees(gamma + PI));
    let zenith = radians_to_degrees(z);

    SolarPosition::new(azimuth, zenith)
}

/// Calculate t parameter (days since 2000-01-01 12:00:00 TT)
fn calc_t<Tz: TimeZone>(datetime: &DateTime<Tz>) -> f64 {
    // Convert to UTC for proper astronomical calculations
    let utc_datetime = datetime.with_timezone(&chrono::Utc);
    let mut m = i32::try_from(utc_datetime.month()).expect("month should fit in i32");
    let mut y = utc_datetime.year();
    let d = i32::try_from(utc_datetime.day()).expect("day should fit in i32");
    let h = f64::from(utc_datetime.hour())
        + f64::from(utc_datetime.minute()) / 60.0
        + f64::from(utc_datetime.second()) / 3600.0;

    if m <= 2 {
        m += 12;
        y -= 1;
    }

    f64::from((365.25 * f64::from(y - 2000)) as i32)
        + f64::from((30.6001 * f64::from(m + 1)) as i32)
        - f64::from((0.01 * f64::from(y)) as i32)
        + f64::from(d)
        + 0.0416667 * h
        - 21958.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{DateTime, FixedOffset};

    #[test]
    fn test_grena3_basic_functionality() {
        let datetime = "2023-06-21T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();

        let result = solar_position(datetime, 37.7749, -122.4194, 69.0, None);

        assert!(result.is_ok());
        let position = result.unwrap();
        assert!(position.azimuth() >= 0.0 && position.azimuth() <= 360.0);
        assert!(position.zenith_angle() >= 0.0 && position.zenith_angle() <= 180.0);
    }

    #[test]
    fn test_grena3_with_refraction() {
        let datetime = "2023-06-21T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();

        let result = solar_position(
            datetime,
            37.7749,
            -122.4194,
            69.0,
            Some(RefractionCorrection::new(1013.25, 15.0).unwrap()),
        );

        assert!(result.is_ok());
        let position = result.unwrap();
        assert!(position.azimuth() >= 0.0 && position.azimuth() <= 360.0);
        assert!(position.zenith_angle() >= 0.0 && position.zenith_angle() <= 180.0);
    }

    #[test]
    fn test_grena3_coordinate_validation() {
        let datetime = "2023-06-21T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();

        // Invalid latitude
        assert!(solar_position(datetime, 95.0, 0.0, 0.0, None).is_err());

        // Invalid longitude
        assert!(solar_position(datetime, 0.0, 185.0, 0.0, None).is_err());
    }

    #[test]
    fn test_calc_t() {
        // Test with a known date
        let datetime = "2023-06-21T12:00:00-07:00"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let t = calc_t(&datetime);

        // The Grena3 algorithm uses a specific reference point that may result in negative values
        // This is correct behavior - just ensure the calculation is consistent
        assert!(t.is_finite(), "t should be finite");

        // Test that the calculation is stable
        let t2 = calc_t(&datetime);
        assert!(
            (t - t2).abs() < f64::EPSILON,
            "calc_t should be deterministic"
        );

        // Test that different dates give different results
        let datetime2 = "2023-06-22T12:00:00Z"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let t3 = calc_t(&datetime2);
        assert!(
            (t - t3).abs() > 0.5,
            "Different dates should give different t values"
        );
    }
}
