//! SPA algorithm implementation.
//!
//! High-accuracy solar positioning based on the NREL algorithm by Reda & Andreas (2003).
//! Accuracy: ±0.0003° for years -2000 to 6000.
//!
//! Reference: Reda, I.; Andreas, A. (2003). Solar position algorithm for solar radiation applications.
//! Solar Energy, 76(5), 577-589. DOI: <http://dx.doi.org/10.1016/j.solener.2003.12.003>

#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::unreadable_literal)]

use crate::error::check_coordinates;
use crate::math::{
    acos, asin, atan, atan2, cos, degrees_to_radians, floor, mul_add, normalize_degrees_0_to_360,
    polynomial, powi, radians_to_degrees, sin, tan,
};
use crate::time::JulianDate;
#[cfg(feature = "chrono")]
use crate::Horizon;
use crate::{RefractionCorrection, Result, SolarPosition};

pub mod coefficients;
use coefficients::{
    NUTATION_COEFFS, OBLIQUITY_COEFFS, TERMS_B, TERMS_L, TERMS_PE, TERMS_R, TERMS_Y,
};

#[cfg(feature = "chrono")]
use chrono::{DateTime, TimeZone, Utc};

/// Standard sunrise/sunset elevation angle (accounts for refraction and sun's radius).
const SUNRISE_SUNSET_ANGLE: f64 = -0.83337;

/// Aberration constant in arcseconds.
const ABERRATION_CONSTANT: f64 = -20.4898;

/// Earth flattening factor (WGS84).
const EARTH_FLATTENING_FACTOR: f64 = 0.99664719;

/// Earth radius in meters (WGS84).
const EARTH_RADIUS_METERS: f64 = 6378140.0;

/// Seconds per hour conversion factor.
const SECONDS_PER_HOUR: f64 = 3600.0;

/// Calculate solar position using the SPA algorithm.
///
/// # Arguments
/// * `datetime` - Date and time with timezone
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `elevation` - Observer elevation in meters above sea level
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `refraction` - Optional atmospheric refraction correction
///
/// # Returns
/// Solar position or error
///
/// # Errors
/// Returns error for invalid coordinates (latitude outside ±90°, longitude outside ±180°)
///
/// # Example
/// ```rust
/// use solar_positioning::{spa, RefractionCorrection};
/// use chrono::{DateTime, FixedOffset};
///
/// let datetime = "2023-06-21T12:00:00-07:00".parse::<DateTime<FixedOffset>>().unwrap();
///
/// // With atmospheric refraction correction
/// let position = spa::solar_position(
///     datetime,
///     37.7749,     // San Francisco latitude
///     -122.4194,   // San Francisco longitude
///     0.0,         // elevation (meters)
///     69.0,        // deltaT (seconds)
///     Some(RefractionCorrection::standard()),
/// ).unwrap();
///
/// // Without refraction correction
/// let position_no_refraction = spa::solar_position(
///     datetime,
///     37.7749,
///     -122.4194,
///     0.0,
///     69.0,
///     None,
/// ).unwrap();
///
/// println!("Azimuth: {:.3}°", position.azimuth());
/// println!("Elevation: {:.3}°", position.elevation_angle());
/// ```
#[cfg(feature = "chrono")]
#[allow(clippy::needless_pass_by_value)]
pub fn solar_position<Tz: TimeZone>(
    datetime: DateTime<Tz>,
    latitude: f64,
    longitude: f64,
    elevation: f64,
    delta_t: f64,
    refraction: Option<RefractionCorrection>,
) -> Result<SolarPosition> {
    let jd = JulianDate::from_datetime(&datetime, delta_t)?;
    solar_position_from_julian(jd, latitude, longitude, elevation, refraction)
}

/// Calculate solar position from a Julian date.
///
/// Core implementation for `no_std` compatibility (no chrono dependency).
///
/// # Arguments
/// * `jd` - Julian date with `delta_t`
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `elevation` - Observer elevation in meters above sea level
/// * `refraction` - Optional atmospheric refraction correction
///
/// # Returns
/// Solar position or error
///
/// # Errors
/// Returns error for invalid coordinates
///
/// # Example
/// ```rust
/// use solar_positioning::{spa, time::JulianDate, RefractionCorrection};
///
/// // Julian date for 2023-06-21 12:00:00 UTC with ΔT=69s
/// let jd = JulianDate::from_utc(2023, 6, 21, 12, 0, 0.0, 69.0).unwrap();
///
/// let position = spa::solar_position_from_julian(
///     jd,
///     37.7749,     // San Francisco latitude
///     -122.4194,   // San Francisco longitude
///     0.0,         // elevation (meters)
///     Some(RefractionCorrection::standard()),
/// ).unwrap();
///
/// println!("Azimuth: {:.3}°", position.azimuth());
/// println!("Elevation: {:.3}°", position.elevation_angle());
/// ```
pub fn solar_position_from_julian(
    jd: JulianDate,
    latitude: f64,
    longitude: f64,
    elevation: f64,
    refraction: Option<RefractionCorrection>,
) -> Result<SolarPosition> {
    let time_dependent = spa_time_dependent_from_julian(jd)?;
    spa_with_time_dependent_parts(latitude, longitude, elevation, refraction, &time_dependent)
}

/// Time-dependent intermediate values from SPA calculation (steps 1-11).
///
/// Pre-computed astronomical values independent of observer location.
/// Use with [`spa_with_time_dependent_parts`] for efficient coordinate sweeps.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SpaTimeDependent {
    /// Geocentric longitude (degrees)
    pub(crate) theta_degrees: f64,
    /// Geocentric latitude (degrees)
    pub(crate) beta_degrees: f64,
    /// Earth radius vector (AU)
    pub(crate) r: f64,
    /// Nutation in longitude (degrees)
    pub(crate) delta_psi: f64,
    /// True obliquity of ecliptic (degrees)
    pub(crate) epsilon_degrees: f64,
    /// Apparent sun longitude (degrees)
    pub(crate) lambda_degrees: f64,
    /// Apparent sidereal time at Greenwich (degrees)
    pub(crate) nu_degrees: f64,
    /// Geocentric sun right ascension (degrees)
    pub(crate) alpha_degrees: f64,
    /// Geocentric sun declination (degrees)
    pub(crate) delta_degrees: f64,
}

#[derive(Debug, Clone, Copy)]
struct DeltaPsiEpsilon {
    delta_psi: f64,
    delta_epsilon: f64,
}

/// Calculate L, B, R terms from the coefficient tables.
fn calculate_lbr_terms(jme: f64, term_coeffs: &[&[&[f64; 3]]]) -> [f64; 6] {
    // We know from coefficients that we have exactly 6 terms for L, 2 for B, and 5 for R
    // Use a fixed-size array to avoid heap allocation
    let mut lbr_terms = [0.0; 6];

    for (i, term_set) in term_coeffs.iter().enumerate().take(6) {
        let mut lbr_sum = 0.0;
        for term in *term_set {
            lbr_sum += term[0] * cos(mul_add(term[2], jme, term[1]));
        }
        lbr_terms[i] = lbr_sum;
    }

    lbr_terms
}

/// Calculate L, B, or R polynomial from the terms.
fn calculate_lbr_polynomial(jme: f64, terms: &[f64], num_terms: usize) -> f64 {
    polynomial(&terms[..num_terms], jme) / 1e8
}

/// Calculate normalized degrees from LBR polynomial
fn lbr_to_normalized_degrees(jme: f64, terms: &[f64], num_terms: usize) -> f64 {
    normalize_degrees_0_to_360(radians_to_degrees(calculate_lbr_polynomial(
        jme, terms, num_terms,
    )))
}

/// Calculate nutation terms (X values).
fn calculate_nutation_terms(jce: f64) -> [f64; 5] {
    // Use fixed-size array to avoid heap allocation
    // NUTATION_COEFFS always has exactly 5 elements
    [
        polynomial(NUTATION_COEFFS[0], jce),
        polynomial(NUTATION_COEFFS[1], jce),
        polynomial(NUTATION_COEFFS[2], jce),
        polynomial(NUTATION_COEFFS[3], jce),
        polynomial(NUTATION_COEFFS[4], jce),
    ]
}

/// Calculate nutation in longitude and obliquity.
fn calculate_delta_psi_epsilon(jce: f64, x: &[f64]) -> DeltaPsiEpsilon {
    let mut delta_psi = 0.0;
    let mut delta_epsilon = 0.0;

    for (i, pe_term) in TERMS_PE.iter().enumerate() {
        let xj_yterm_sum = degrees_to_radians(calculate_xj_yterm_sum(i, x));

        // Use Math.fma equivalent: a * b + c
        let delta_psi_contrib = mul_add(pe_term[1], jce, pe_term[0]) * sin(xj_yterm_sum);
        let delta_epsilon_contrib = mul_add(pe_term[3], jce, pe_term[2]) * cos(xj_yterm_sum);

        delta_psi += delta_psi_contrib;
        delta_epsilon += delta_epsilon_contrib;
    }

    DeltaPsiEpsilon {
        delta_psi: delta_psi / 36_000_000.0,
        delta_epsilon: delta_epsilon / 36_000_000.0,
    }
}

/// Calculate sum of X[j] * Y[i][j] for nutation.
fn calculate_xj_yterm_sum(i: usize, x: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (j, &x_val) in x.iter().enumerate() {
        sum += x_val * f64::from(TERMS_Y[i][j]);
    }
    sum
}

/// Calculate true obliquity of the ecliptic.
fn calculate_true_obliquity_of_ecliptic(jd: &JulianDate, delta_epsilon: f64) -> f64 {
    let epsilon0 = polynomial(OBLIQUITY_COEFFS, jd.julian_ephemeris_millennium() / 10.0);
    epsilon0 / 3600.0 + delta_epsilon
}

/// Calculate apparent sidereal time at Greenwich.
fn calculate_apparent_sidereal_time_at_greenwich(
    jd: &JulianDate,
    delta_psi: f64,
    epsilon_degrees: f64,
) -> f64 {
    let nu0_degrees = normalize_degrees_0_to_360(mul_add(
        powi(jd.julian_century(), 2),
        0.000387933 - jd.julian_century() / 38710000.0,
        mul_add(
            360.98564736629f64,
            jd.julian_date() - 2451545.0,
            280.46061837,
        ),
    ));

    mul_add(
        delta_psi,
        cos(degrees_to_radians(epsilon_degrees)),
        nu0_degrees,
    )
}

/// Calculate geocentric sun right ascension.
fn calculate_geocentric_sun_right_ascension(
    beta_rad: f64,
    epsilon_rad: f64,
    lambda_rad: f64,
) -> f64 {
    let alpha = atan2(
        mul_add(
            sin(lambda_rad),
            cos(epsilon_rad),
            -(tan(beta_rad) * sin(epsilon_rad)),
        ),
        cos(lambda_rad),
    );
    normalize_degrees_0_to_360(radians_to_degrees(alpha))
}

/// Calculate geocentric sun declination.
fn calculate_geocentric_sun_declination(beta_rad: f64, epsilon_rad: f64, lambda_rad: f64) -> f64 {
    asin(mul_add(
        sin(beta_rad),
        cos(epsilon_rad),
        cos(beta_rad) * sin(epsilon_rad) * sin(lambda_rad),
    ))
}

/// Calculate sunrise/sunset times without chrono dependency.
///
/// Returns times as hours since midnight UTC (0.0 to 24.0+) for the given date.
/// Hours can extend beyond 24.0 (next day) or be negative (previous day).
///
/// This follows the NREL SPA algorithm (Reda & Andreas 2003, Appendix A.2).
///
/// # Arguments
/// * `year` - Year (can be negative for BCE)
/// * `month` - Month (1-12)
/// * `day` - Day of month (1-31)
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `elevation_angle` - Sun elevation angle for sunrise/sunset in degrees (typically -0.833°)
///
/// # Returns
/// `SunriseResult<HoursUtc>` with times as hours since midnight UTC
///
/// # Errors
/// Returns error for invalid date components or coordinates
///
/// # Example
/// ```
/// use solar_positioning::{spa, HoursUtc};
///
/// let result = spa::sunrise_sunset_utc(
///     2023, 6, 21,   // June 21, 2023
///     37.7749,       // San Francisco latitude
///     -122.4194,     // San Francisco longitude
///     69.0,          // deltaT (seconds)
///     -0.833         // standard sunrise/sunset angle
/// ).unwrap();
///
/// if let solar_positioning::SunriseResult::RegularDay { sunrise, transit, sunset } = result {
///     println!("Sunrise: {:.2} hours UTC", sunrise.hours());
///     println!("Transit: {:.2} hours UTC", transit.hours());
///     println!("Sunset: {:.2} hours UTC", sunset.hours());
/// }
/// ```
pub fn sunrise_sunset_utc(
    year: i32,
    month: u32,
    day: u32,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    elevation_angle: f64,
) -> Result<crate::SunriseResult<crate::HoursUtc>> {
    check_coordinates(latitude, longitude)?;

    // Create Julian date for midnight UTC (0 UT) of the given date
    let jd_midnight = JulianDate::from_utc(year, month, day, 0, 0, 0.0, 0.0)?;

    // Calculate sunrise/sunset using core algorithm
    calculate_sunrise_sunset_core(jd_midnight, latitude, longitude, delta_t, elevation_angle)
}

/// Calculate sunrise, solar transit, and sunset times for a specific horizon type.
///
/// This is a convenience function that uses predefined elevation angles for common
/// sunrise/twilight calculations without requiring the chrono library.
///
/// # Arguments
/// * `year` - Year (e.g., 2023)
/// * `month` - Month (1-12)
/// * `day` - Day of month (1-31)
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `horizon` - Horizon type (sunrise/sunset, civil twilight, etc.)
///
/// # Returns
/// `SunriseResult<HoursUtc>` with times as hours since midnight UTC
///
/// # Errors
/// Returns error for invalid coordinates, dates, or computation failures.
///
/// # Example
/// ```rust
/// use solar_positioning::{spa, Horizon};
///
/// // Standard sunrise/sunset
/// let result = spa::sunrise_sunset_utc_for_horizon(
///     2023, 6, 21,
///     37.7749,   // San Francisco latitude
///     -122.4194, // San Francisco longitude
///     69.0,      // deltaT (seconds)
///     Horizon::SunriseSunset
/// ).unwrap();
///
/// // Civil twilight
/// let twilight = spa::sunrise_sunset_utc_for_horizon(
///     2023, 6, 21,
///     37.7749, -122.4194, 69.0,
///     Horizon::CivilTwilight
/// ).unwrap();
/// ```
pub fn sunrise_sunset_utc_for_horizon(
    year: i32,
    month: u32,
    day: u32,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    horizon: crate::Horizon,
) -> Result<crate::SunriseResult<crate::HoursUtc>> {
    sunrise_sunset_utc(
        year,
        month,
        day,
        latitude,
        longitude,
        delta_t,
        horizon.elevation_angle(),
    )
}

/// Calculate sunrise, solar transit, and sunset times using the SPA algorithm.
///
/// This follows the NREL SPA algorithm (Reda & Andreas 2003) for calculating
/// sunrise, transit (solar noon), and sunset times with high accuracy.
///
/// # Arguments
/// * `date` - Date for calculations
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `elevation_angle` - Sun elevation angle for sunrise/sunset in degrees (typically -0.833°)
///
/// # Returns
/// `SunriseResult` variant indicating regular day, polar day, or polar night
///
/// # Errors
/// Returns error for invalid coordinates (latitude outside ±90°, longitude outside ±180°)
///
/// # Panics
/// May panic if date cannot be converted to start of day (unlikely with valid dates).
///
/// # Example
/// ```rust
/// use solar_positioning::spa;
/// use chrono::{DateTime, FixedOffset, NaiveDate, TimeZone};
///
/// let date = FixedOffset::east_opt(-7 * 3600).unwrap() // Pacific Time (UTC-7)
///     .from_local_datetime(&NaiveDate::from_ymd_opt(2023, 6, 21).unwrap()
///         .and_hms_opt(0, 0, 0).unwrap()).unwrap();
/// let result = spa::sunrise_sunset(
///     date,
///     37.7749,   // San Francisco latitude
///     -122.4194, // San Francisco longitude
///     69.0,      // deltaT (seconds)
///     -0.833     // standard sunrise/sunset angle
/// ).unwrap();
#[cfg(feature = "chrono")]
#[allow(clippy::needless_pass_by_value)]
pub fn sunrise_sunset<Tz: TimeZone>(
    date: DateTime<Tz>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    elevation_angle: f64,
) -> Result<crate::SunriseResult<DateTime<Tz>>> {
    check_coordinates(latitude, longitude)?;

    let day_start = truncate_to_day_start(&date);

    // Implement the complete algorithm in one place for accuracy
    calculate_sunrise_sunset_spa_algorithm(day_start, latitude, longitude, delta_t, elevation_angle)
}

/// Core sunrise/sunset calculation that returns times as fractions of day.
///
/// This is the shared implementation used by both chrono and non-chrono APIs.
#[allow(clippy::unnecessary_wraps)]
fn calculate_sunrise_sunset_core(
    jd_midnight: JulianDate,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    elevation_angle: f64,
) -> Result<crate::SunriseResult<crate::HoursUtc>> {
    // A.2.1. Calculate the apparent sidereal time at Greenwich at 0 UT
    let jce_day = jd_midnight.julian_ephemeris_century();
    let x_terms = calculate_nutation_terms(jce_day);
    let delta_psi_epsilon = calculate_delta_psi_epsilon(jce_day, &x_terms);
    let epsilon_degrees =
        calculate_true_obliquity_of_ecliptic(&jd_midnight, delta_psi_epsilon.delta_epsilon);
    let nu_degrees = calculate_apparent_sidereal_time_at_greenwich(
        &jd_midnight,
        delta_psi_epsilon.delta_psi,
        epsilon_degrees,
    );

    // A.2.2. Calculate alpha/delta for day before, same day, next day
    let mut alpha_deltas = [AlphaDelta {
        alpha: 0.0,
        delta: 0.0,
    }; 3];
    for (i, alpha_delta) in alpha_deltas.iter_mut().enumerate() {
        let current_jd = jd_midnight.add_days((i as f64) - 1.0);
        let current_jme = current_jd.julian_ephemeris_millennium();
        let ad = calculate_alpha_delta(current_jme, delta_psi_epsilon.delta_psi, epsilon_degrees);
        *alpha_delta = ad;
    }

    // Calculate initial transit time and check for polar conditions
    let m0 = (alpha_deltas[1].alpha - longitude - nu_degrees) / 360.0;

    // Check for polar conditions but don't return early - need to calculate accurate transit
    let polar_type = check_polar_conditions_type(latitude, elevation_angle, alpha_deltas[1].delta);

    // Calculate approximate times and apply corrections (needed for accurate transit time)
    let (m_values, _h0_degrees) =
        calculate_approximate_times(m0, latitude, elevation_angle, alpha_deltas[1].delta);

    // Apply final corrections to get accurate times (as fractions of day)
    let (t_frac, r_frac, s_frac) = calculate_final_time_fractions(
        m_values,
        nu_degrees,
        delta_t,
        latitude,
        longitude,
        elevation_angle,
        alpha_deltas,
    );

    // Convert fractions to hours (0-24+, can be negative or > 24)
    let transit_hours = crate::HoursUtc::from_hours(t_frac * 24.0);
    let sunrise_hours = crate::HoursUtc::from_hours(r_frac * 24.0);
    let sunset_hours = crate::HoursUtc::from_hours(s_frac * 24.0);

    // Return appropriate result type based on polar conditions
    match polar_type {
        Some(PolarType::AllDay) => Ok(crate::SunriseResult::AllDay {
            transit: transit_hours,
        }),
        Some(PolarType::AllNight) => Ok(crate::SunriseResult::AllNight {
            transit: transit_hours,
        }),
        None => Ok(crate::SunriseResult::RegularDay {
            sunrise: sunrise_hours,
            transit: transit_hours,
            sunset: sunset_hours,
        }),
    }
}

/// Calculate sunrise/sunset times using NREL SPA algorithm Appendix A.2
#[cfg(feature = "chrono")]
fn calculate_sunrise_sunset_spa_algorithm<Tz: TimeZone>(
    day: DateTime<Tz>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    elevation_angle: f64,
) -> Result<crate::SunriseResult<DateTime<Tz>>> {
    let day_start = truncate_to_day_start(&day);

    // A.2.1. Calculate the apparent sidereal time at Greenwich at 0 UT
    let (nu_degrees, delta_psi_epsilon, epsilon_degrees) =
        calculate_sidereal_time_and_nutation(day_start.clone());

    // A.2.2. Calculate alpha/delta for day before, same day, next day
    let alpha_deltas =
        calculate_alpha_deltas_for_three_days(day_start, delta_psi_epsilon, epsilon_degrees)?;

    // Calculate initial transit time and check for polar conditions
    let m0 = (alpha_deltas[1].alpha - longitude - nu_degrees) / 360.0;

    // Check for polar conditions but don't return early - need to calculate accurate transit
    let polar_type = check_polar_conditions_type(latitude, elevation_angle, alpha_deltas[1].delta);

    // Calculate approximate times and apply corrections (needed for accurate transit time)
    let (m_values, _h0_degrees) =
        calculate_approximate_times(m0, latitude, elevation_angle, alpha_deltas[1].delta);

    // Apply final corrections to get accurate times
    let final_result = calculate_final_times(FinalTimeParams {
        day,
        m_values,
        nu_degrees,
        delta_t,
        latitude,
        longitude,
        elevation_angle,
        alpha_deltas,
    });

    // Convert to appropriate polar result type if needed
    match polar_type {
        Some(PolarType::AllDay) => {
            if let crate::SunriseResult::RegularDay { transit, .. } = final_result {
                Ok(crate::SunriseResult::AllDay { transit })
            } else {
                Ok(final_result)
            }
        }
        Some(PolarType::AllNight) => {
            if let crate::SunriseResult::RegularDay { transit, .. } = final_result {
                Ok(crate::SunriseResult::AllNight { transit })
            } else {
                Ok(final_result)
            }
        }
        None => Ok(final_result),
    }
}
#[cfg(feature = "chrono")]
#[allow(clippy::needless_pass_by_value)]
fn calculate_sidereal_time_and_nutation<Tz: TimeZone>(
    day_start: DateTime<Tz>,
) -> (f64, DeltaPsiEpsilon, f64) {
    let jd = JulianDate::from_datetime(&day_start, 0.0).unwrap();
    let jce_day = jd.julian_ephemeris_century();
    let x_terms = calculate_nutation_terms(jce_day);
    let delta_psi_epsilon = calculate_delta_psi_epsilon(jce_day, &x_terms);
    let epsilon_degrees =
        calculate_true_obliquity_of_ecliptic(&jd, delta_psi_epsilon.delta_epsilon);
    let nu_degrees = calculate_apparent_sidereal_time_at_greenwich(
        &jd,
        delta_psi_epsilon.delta_psi,
        epsilon_degrees,
    );
    (nu_degrees, delta_psi_epsilon, epsilon_degrees)
}

#[cfg(feature = "chrono")]
#[allow(clippy::needless_pass_by_value)]
fn calculate_alpha_deltas_for_three_days<Tz: TimeZone>(
    day_start: DateTime<Tz>,
    delta_psi_epsilon: DeltaPsiEpsilon,
    epsilon_degrees: f64,
) -> Result<[AlphaDelta; 3]> {
    let mut alpha_deltas = [AlphaDelta {
        alpha: 0.0,
        delta: 0.0,
    }; 3];
    for (i, alpha_delta) in alpha_deltas.iter_mut().enumerate() {
        let current_jd = JulianDate::from_datetime(&day_start, 0.0)?.add_days((i as f64) - 1.0);
        let current_jme = current_jd.julian_ephemeris_millennium();
        let ad = calculate_alpha_delta(current_jme, delta_psi_epsilon.delta_psi, epsilon_degrees);
        *alpha_delta = ad;
    }
    Ok(alpha_deltas)
}

// ============================================================================
// Sunrise/sunset helper functions below
// Core functions work without chrono, chrono-specific wrappers separate
// ============================================================================

/// Enum for polar condition types
#[derive(Debug, Clone, Copy)]
enum PolarType {
    AllDay,
    AllNight,
}

/// Check for polar day/night conditions and return the type
fn check_polar_conditions_type(
    latitude: f64,
    elevation_angle: f64,
    delta1: f64,
) -> Option<PolarType> {
    let phi = degrees_to_radians(latitude);
    let elevation_rad = degrees_to_radians(elevation_angle);
    let delta1_rad = degrees_to_radians(delta1);

    let acos_arg =
        mul_add(sin(phi), -sin(delta1_rad), sin(elevation_rad)) / (cos(phi) * cos(delta1_rad));

    if acos_arg < -1.0 {
        Some(PolarType::AllDay)
    } else if acos_arg > 1.0 {
        Some(PolarType::AllNight)
    } else {
        None
    }
}

/// Check for polar day/night conditions and construct `SunriseResult`
#[cfg(feature = "chrono")]
fn check_polar_conditions<Tz: TimeZone>(
    day: DateTime<Tz>,
    m0: f64,
    latitude: f64,
    elevation_angle: f64,
    delta1: f64,
) -> Option<crate::SunriseResult<DateTime<Tz>>> {
    check_polar_conditions_type(latitude, elevation_angle, delta1).map(|polar_type| {
        let transit = add_fraction_of_day(day, m0);
        match polar_type {
            PolarType::AllDay => crate::SunriseResult::AllDay { transit },
            PolarType::AllNight => crate::SunriseResult::AllNight { transit },
        }
    })
}

/// A.2.5-6. Calculate approximate times for transit, sunrise, sunset
fn calculate_approximate_times(
    m0: f64,
    latitude: f64,
    elevation_angle: f64,
    delta1: f64,
) -> ([f64; 3], f64) {
    let phi = degrees_to_radians(latitude);
    let delta1_rad = degrees_to_radians(delta1);
    let elevation_rad = degrees_to_radians(elevation_angle);

    let acos_arg =
        mul_add(sin(phi), -sin(delta1_rad), sin(elevation_rad)) / (cos(phi) * cos(delta1_rad));
    let h0 = acos(acos_arg);
    let h0_degrees = radians_to_degrees(h0).min(180.0);

    let mut m = [0.0; 3];
    m[0] = normalize_to_unit_range(m0);
    m[1] = normalize_to_unit_range(m0 - h0_degrees / 360.0);
    m[2] = normalize_to_unit_range(m0 + h0_degrees / 360.0);

    (m, h0_degrees)
}

/// A.2.8-15. Calculate final accurate time fractions using corrections
/// Returns (`transit_frac`, `sunrise_frac`, `sunset_frac`) as fractions of day
fn calculate_final_time_fractions(
    m_values: [f64; 3],
    nu_degrees: f64,
    delta_t: f64,
    latitude: f64,
    longitude: f64,
    elevation_angle: f64,
    alpha_deltas: [AlphaDelta; 3],
) -> (f64, f64, f64) {
    // A.2.8. Calculate sidereal times
    let mut nu = [0.0; 3];
    for (i, nu_item) in nu.iter_mut().enumerate() {
        *nu_item = mul_add(360.985647f64, m_values[i], nu_degrees);
    }

    // A.2.9. Calculate terms with deltaT correction
    let mut n = [0.0; 3];
    for (i, n_item) in n.iter_mut().enumerate() {
        *n_item = m_values[i] + delta_t / 86400.0;
    }

    // A.2.10. Calculate α'i and δ'i using interpolation
    let alpha_delta_primes = calculate_interpolated_alpha_deltas(&alpha_deltas, &n);

    // A.2.11. Calculate local hour angles
    let mut h_prime = [0.0; 3];
    for i in 0..3 {
        let h_prime_i = nu[i] + longitude - alpha_delta_primes[i].alpha;
        h_prime[i] = limit_h_prime(h_prime_i);
    }

    // A.2.12. Calculate sun altitudes
    let phi = degrees_to_radians(latitude);
    let mut h = [0.0; 3];
    for i in 0..3 {
        let delta_prime_rad = degrees_to_radians(alpha_delta_primes[i].delta);
        h[i] = radians_to_degrees(asin(mul_add(
            sin(phi),
            sin(delta_prime_rad),
            cos(phi) * cos(delta_prime_rad) * cos(degrees_to_radians(h_prime[i])),
        )));
    }

    // A.2.13-15. Calculate final times as fractions
    let t = m_values[0] - h_prime[0] / 360.0;
    let r = m_values[1]
        + (h[1] - elevation_angle)
            / (360.0
                * cos(degrees_to_radians(alpha_delta_primes[1].delta))
                * cos(phi)
                * sin(degrees_to_radians(h_prime[1])));
    let s = m_values[2]
        + (h[2] - elevation_angle)
            / (360.0
                * cos(degrees_to_radians(alpha_delta_primes[2].delta))
                * cos(phi)
                * sin(degrees_to_radians(h_prime[2])));

    (t, r, s)
}

/// A.2.8-15. Calculate final accurate times using corrections
#[cfg(feature = "chrono")]
fn calculate_final_times<Tz: TimeZone>(
    params: FinalTimeParams<Tz>,
) -> crate::SunriseResult<DateTime<Tz>> {
    // A.2.8. Calculate sidereal times
    let mut nu = [0.0; 3];
    for (i, nu_item) in nu.iter_mut().enumerate() {
        *nu_item = mul_add(360.985647f64, params.m_values[i], params.nu_degrees);
    }

    // A.2.9. Calculate terms with deltaT correction
    let mut n = [0.0; 3];
    for (i, n_item) in n.iter_mut().enumerate() {
        *n_item = params.m_values[i] + params.delta_t / 86400.0;
    }

    // A.2.10. Calculate α'i and δ'i using interpolation
    let alpha_delta_primes = calculate_interpolated_alpha_deltas(&params.alpha_deltas, &n);

    // A.2.11. Calculate local hour angles
    let mut h_prime = [0.0; 3];
    for i in 0..3 {
        let h_prime_i = nu[i] + params.longitude - alpha_delta_primes[i].alpha;
        h_prime[i] = limit_h_prime(h_prime_i);
    }

    // A.2.12. Calculate sun altitudes
    let phi = degrees_to_radians(params.latitude);
    let mut h = [0.0; 3];
    for i in 0..3 {
        let delta_prime_rad = degrees_to_radians(alpha_delta_primes[i].delta);
        h[i] = radians_to_degrees(asin(mul_add(
            sin(phi),
            sin(delta_prime_rad),
            cos(phi) * cos(delta_prime_rad) * cos(degrees_to_radians(h_prime[i])),
        )));
    }

    // A.2.13-15. Calculate final times
    let t = params.m_values[0] - h_prime[0] / 360.0;
    let r = params.m_values[1]
        + (h[1] - params.elevation_angle)
            / (360.0
                * cos(degrees_to_radians(alpha_delta_primes[1].delta))
                * cos(phi)
                * sin(degrees_to_radians(h_prime[1])));
    let s = params.m_values[2]
        + (h[2] - params.elevation_angle)
            / (360.0
                * cos(degrees_to_radians(alpha_delta_primes[2].delta))
                * cos(phi)
                * sin(degrees_to_radians(h_prime[2])));

    let sunrise = add_fraction_of_day(params.day.clone(), r);
    let transit = add_fraction_of_day(params.day.clone(), t);
    let sunset = add_fraction_of_day(params.day, s);

    crate::SunriseResult::RegularDay {
        sunrise,
        transit,
        sunset,
    }
}

/// A.2.10. Calculate interpolated alpha/delta values
fn calculate_interpolated_alpha_deltas(
    alpha_deltas: &[AlphaDelta; 3],
    n: &[f64; 3],
) -> [AlphaDelta; 3] {
    let a = limit_if_necessary(alpha_deltas[1].alpha - alpha_deltas[0].alpha);
    let a_prime = limit_if_necessary(alpha_deltas[1].delta - alpha_deltas[0].delta);

    let b = limit_if_necessary(alpha_deltas[2].alpha - alpha_deltas[1].alpha);
    let b_prime = limit_if_necessary(alpha_deltas[2].delta - alpha_deltas[1].delta);

    let c = b - a;
    let c_prime = b_prime - a_prime;

    let mut alpha_delta_primes = [AlphaDelta {
        alpha: 0.0,
        delta: 0.0,
    }; 3];
    for i in 0..3 {
        alpha_delta_primes[i].alpha =
            alpha_deltas[1].alpha + (n[i] * (mul_add(c, n[i], a + b))) / 2.0;
        alpha_delta_primes[i].delta =
            alpha_deltas[1].delta + (n[i] * (mul_add(c_prime, n[i], a_prime + b_prime))) / 2.0;
    }
    alpha_delta_primes
}

#[derive(Debug, Clone, Copy)]
struct AlphaDelta {
    alpha: f64,
    delta: f64,
}

/// Parameters for calculating final sunrise/sunset times
#[derive(Debug, Clone)]
#[cfg(feature = "chrono")]
struct FinalTimeParams<Tz: TimeZone> {
    day: DateTime<Tz>,
    m_values: [f64; 3],
    nu_degrees: f64,
    delta_t: f64,
    latitude: f64,
    longitude: f64,
    elevation_angle: f64,
    alpha_deltas: [AlphaDelta; 3],
}

/// Calculate sunrise, solar transit, and sunset times for a specific horizon type.
///
/// This is a convenience function that uses predefined elevation angles for common
/// sunrise/twilight calculations.
///
/// # Arguments
/// * `date` - Date for calculations
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `horizon` - Horizon type (sunrise/sunset, civil twilight, etc.)
///
/// # Errors
/// Returns error for invalid coordinates, dates, or computation failures.
///
/// # Panics
/// May panic if date cannot be converted to start of day (unlikely with valid dates).
///
/// # Example
/// ```rust
/// use solar_positioning::{spa, Horizon};
/// use chrono::{FixedOffset, NaiveDate, TimeZone};
///
/// let date = FixedOffset::east_opt(-7 * 3600).unwrap() // Pacific Time (UTC-7)
///     .from_local_datetime(&NaiveDate::from_ymd_opt(2023, 6, 21).unwrap()
///         .and_hms_opt(0, 0, 0).unwrap()).unwrap();
///
/// // Standard sunrise/sunset
/// let sunrise_result = spa::sunrise_sunset_for_horizon(
///     date, 37.7749, -122.4194, 69.0, Horizon::SunriseSunset
/// ).unwrap();
///
/// // Civil twilight
/// let twilight_result = spa::sunrise_sunset_for_horizon(
///     date, 37.7749, -122.4194, 69.0, Horizon::CivilTwilight
/// ).unwrap();
/// ```
#[cfg(feature = "chrono")]
pub fn sunrise_sunset_for_horizon<Tz: TimeZone>(
    date: DateTime<Tz>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    horizon: Horizon,
) -> Result<crate::SunriseResult<DateTime<Tz>>> {
    sunrise_sunset(
        date,
        latitude,
        longitude,
        delta_t,
        horizon.elevation_angle(),
    )
}

/// Calculate alpha (right ascension) and delta (declination) for a given JME using full SPA algorithm
/// Following NREL SPA Algorithm Section 3.2-3.8 for sunrise/sunset calculations
fn calculate_alpha_delta(jme: f64, delta_psi: f64, epsilon_degrees: f64) -> AlphaDelta {
    // Follow Java calculateAlphaDelta exactly

    // 3.2.3. Calculate Earth heliocentric latitude, B
    let b_terms = calculate_lbr_terms(jme, TERMS_B);
    let b_degrees = lbr_to_normalized_degrees(jme, &b_terms, TERMS_B.len());

    // 3.2.4. Calculate Earth radius vector, R
    let r_terms = calculate_lbr_terms(jme, TERMS_R);
    let r = calculate_lbr_polynomial(jme, &r_terms, TERMS_R.len());
    assert!(
        r != 0.0,
        "Earth radius vector is zero - astronomical impossibility"
    );

    // 3.2.2. Calculate Earth heliocentric longitude, L
    let l_terms = calculate_lbr_terms(jme, TERMS_L);
    let l_degrees = lbr_to_normalized_degrees(jme, &l_terms, TERMS_L.len());

    // 3.2.5. Calculate geocentric longitude, theta
    let theta_degrees = normalize_degrees_0_to_360(l_degrees + 180.0);

    // 3.2.6. Calculate geocentric latitude, beta
    let beta_degrees = -b_degrees;
    let beta = degrees_to_radians(beta_degrees);
    let epsilon = degrees_to_radians(epsilon_degrees);

    // 3.5. Calculate aberration correction
    let delta_tau = ABERRATION_CONSTANT / (SECONDS_PER_HOUR * r);

    // 3.6. Calculate the apparent sun longitude
    let lambda_degrees = theta_degrees + delta_psi + delta_tau;
    let lambda = degrees_to_radians(lambda_degrees);

    // 3.8.1-3.8.2. Calculate the geocentric sun right ascension and declination
    let alpha_degrees = calculate_geocentric_sun_right_ascension(beta, epsilon, lambda);
    let delta_degrees =
        radians_to_degrees(calculate_geocentric_sun_declination(beta, epsilon, lambda));

    AlphaDelta {
        alpha: alpha_degrees,
        delta: delta_degrees,
    }
}

/// Normalize value to [0, 1) range using the same logic as the removed `limit_to` function
fn normalize_to_unit_range(val: f64) -> f64 {
    let divided = val;
    let limited = divided - floor(divided);
    if limited < 0.0 {
        limited + 1.0
    } else {
        limited
    }
}

#[cfg(feature = "chrono")]
fn truncate_to_day_start<Tz: TimeZone>(datetime: &DateTime<Tz>) -> DateTime<Tz> {
    let tz = datetime.timezone();
    let utc_midnight = datetime
        .with_timezone(&Utc)
        .date_naive()
        .and_hms_opt(0, 0, 0)
        .expect("midnight is always valid")
        .and_utc();

    tz.from_utc_datetime(&utc_midnight.naive_utc())
}

#[cfg(feature = "chrono")]
#[allow(clippy::needless_pass_by_value)]
fn add_fraction_of_day<Tz: TimeZone>(day: DateTime<Tz>, fraction: f64) -> DateTime<Tz> {
    // Match Java implementation exactly:
    // 1. Truncate to start of day (like Java's truncatedTo(ChronoUnit.DAYS))
    // 2. Calculate milliseconds as int (truncating fractional milliseconds)
    // 3. Add the milliseconds
    const MS_PER_DAY: i32 = 24 * 60 * 60 * 1000; // Use i32 like Java
    let millis_plus = (f64::from(MS_PER_DAY) * fraction) as i32; // Cast to i32 like Java

    let day_start = truncate_to_day_start(&day);

    day_start + chrono::Duration::milliseconds(i64::from(millis_plus))
}

/// Limit to 0..1 if absolute value > 2 (Java limitIfNecessary)
fn limit_if_necessary(val: f64) -> f64 {
    if val.abs() > 2.0 {
        normalize_to_unit_range(val)
    } else {
        val
    }
}

/// Limit H' values according to A.2.11
fn limit_h_prime(h_prime: f64) -> f64 {
    let normalized = h_prime / 360.0;
    let limited = 360.0 * (normalized - floor(normalized));

    if limited < -180.0 {
        limited + 360.0
    } else if limited > 180.0 {
        limited - 360.0
    } else {
        limited
    }
}

/// Calculate sunrise/sunset times for multiple horizons efficiently.
///
/// Returns an iterator that yields `(Horizon, SunriseResult)` pairs. This is more
/// efficient than separate calls as it reuses expensive astronomical calculations.
///
/// # Arguments
/// * `date` - Date for calculations
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `horizons` - Iterator of horizon types to calculate
///
/// # Returns
/// Iterator over `Result<(Horizon, SunriseResult)>`
///
/// # Errors
/// Returns error for invalid coordinates (latitude outside ±90°, longitude outside ±180°)
///
/// # Panics
/// May panic if date cannot be converted to start of day (unlikely with valid dates).
///
/// # Example
/// ```rust
/// use solar_positioning::{spa, Horizon};
/// use chrono::{DateTime, FixedOffset};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let datetime = "2023-06-21T12:00:00-07:00".parse::<DateTime<FixedOffset>>().unwrap();
/// let horizons = [
///     Horizon::SunriseSunset,
///     Horizon::CivilTwilight,
///     Horizon::NauticalTwilight,
/// ];
///
/// let results: Result<Vec<_>, _> = spa::sunrise_sunset_multiple(
///     datetime,
///     37.7749,     // San Francisco latitude
///     -122.4194,   // San Francisco longitude
///     69.0,        // deltaT (seconds)
///     horizons.iter().copied()
/// ).collect();
///
/// for (horizon, result) in results? {
///     println!("{:?}: {:?}", horizon, result);
/// }
/// # Ok(())
/// # }
/// ```
#[cfg(feature = "chrono")]
pub fn sunrise_sunset_multiple<Tz, H>(
    date: DateTime<Tz>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    horizons: H,
) -> impl Iterator<Item = Result<(Horizon, crate::SunriseResult<DateTime<Tz>>)>>
where
    Tz: TimeZone,
    H: IntoIterator<Item = Horizon>,
{
    // Pre-calculate common values once for efficiency
    let precomputed = (|| -> Result<_> {
        check_coordinates(latitude, longitude)?;

        let day_start = truncate_to_day_start(&date);
        let (nu_degrees, delta_psi_epsilon, epsilon_degrees) =
            calculate_sidereal_time_and_nutation(day_start.clone());
        let alpha_deltas =
            calculate_alpha_deltas_for_three_days(day_start, delta_psi_epsilon, epsilon_degrees)?;

        Ok((nu_degrees, alpha_deltas))
    })();

    horizons.into_iter().map(move |horizon| match &precomputed {
        Ok((nu_degrees, alpha_deltas)) => {
            let result = calculate_sunrise_sunset_spa_algorithm_with_precomputed(
                date.clone(),
                latitude,
                longitude,
                delta_t,
                horizon.elevation_angle(),
                *nu_degrees,
                *alpha_deltas,
            );
            Ok((horizon, result))
        }
        Err(e) => Err(e.clone()),
    })
}
#[cfg(feature = "chrono")]
#[allow(clippy::needless_pass_by_value)]
fn calculate_sunrise_sunset_spa_algorithm_with_precomputed<Tz: TimeZone>(
    day: DateTime<Tz>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    elevation_angle: f64,
    nu_degrees: f64,
    alpha_deltas: [AlphaDelta; 3],
) -> crate::SunriseResult<DateTime<Tz>> {
    let day_start = truncate_to_day_start(&day);

    // Calculate initial transit time and check for polar conditions
    let m0 = (alpha_deltas[1].alpha - longitude - nu_degrees) / 360.0;
    if let Some(polar_result) = check_polar_conditions(
        day_start.clone(),
        m0,
        latitude,
        elevation_angle,
        alpha_deltas[1].delta,
    ) {
        return polar_result;
    }

    // Calculate approximate times and apply corrections
    let (m_values, _h0_degrees) =
        calculate_approximate_times(m0, latitude, elevation_angle, alpha_deltas[1].delta);

    // Apply final corrections to get accurate times
    calculate_final_times(FinalTimeParams {
        day: day_start,
        m_values,
        nu_degrees,
        delta_t,
        latitude,
        longitude,
        elevation_angle,
        alpha_deltas,
    })
}

/// Extract expensive time-dependent parts of SPA calculation (steps 1-11).
///
/// This function calculates the expensive astronomical quantities that are independent
/// of observer location. Typically used for coordinate sweeps (many locations at fixed
/// time).
///
/// # Arguments
/// * `datetime` - Date and time with timezone
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
///
/// # Returns
/// Pre-computed time-dependent values for SPA calculations
///
/// # Performance
///
/// Use this with [`spa_with_time_dependent_parts`] for coordinate sweeps:
/// ```rust
/// use solar_positioning::spa;
/// use chrono::{DateTime, Utc};
///
/// let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
/// let shared_parts = spa::spa_time_dependent_parts(datetime, 69.0)?;
///
/// for lat in -60..=60 {
///     for lon in -180..=179 {
///         let pos = spa::spa_with_time_dependent_parts(
///             lat as f64, lon as f64, 0.0, None, &shared_parts
///         )?;
///     }
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
/// Returns error if Julian date calculation fails for the provided datetime
///
/// # Panics
#[cfg(feature = "chrono")]
#[allow(clippy::needless_pass_by_value)]
pub fn spa_time_dependent_parts<Tz: TimeZone>(
    datetime: DateTime<Tz>,
    delta_t: f64,
) -> Result<SpaTimeDependent> {
    let jd = JulianDate::from_datetime(&datetime, delta_t)?;
    spa_time_dependent_from_julian(jd)
}

/// Calculate time-dependent parts of SPA from a Julian date.
///
/// Core implementation for `no_std` compatibility.
///
/// # Errors
/// Returns error if Julian date is invalid.
///
/// # Panics
/// Panics if Earth radius vector is zero (astronomical impossibility).
pub fn spa_time_dependent_from_julian(jd: JulianDate) -> Result<SpaTimeDependent> {
    let jme = jd.julian_ephemeris_millennium();
    let jce = jd.julian_ephemeris_century();

    // 3.2.2. Calculate the Earth heliocentric longitude, L (in degrees)
    let l_terms = calculate_lbr_terms(jme, TERMS_L);
    let l_degrees = lbr_to_normalized_degrees(jme, &l_terms, TERMS_L.len());

    // 3.2.3. Calculate the Earth heliocentric latitude, B (in degrees)
    let b_terms = calculate_lbr_terms(jme, TERMS_B);
    let b_degrees = lbr_to_normalized_degrees(jme, &b_terms, TERMS_B.len());

    // 3.2.4. Calculate the Earth radius vector, R (in Astronomical Units, AU)
    let r_terms = calculate_lbr_terms(jme, TERMS_R);
    let r = calculate_lbr_polynomial(jme, &r_terms, TERMS_R.len());

    // Earth's radius vector should never be zero (would mean Earth at center of Sun)
    assert!(
        r != 0.0,
        "Earth radius vector is zero - astronomical impossibility"
    );

    // 3.2.5. Calculate the geocentric longitude, theta (in degrees)
    let theta_degrees = normalize_degrees_0_to_360(l_degrees + 180.0);
    // 3.2.6. Calculate the geocentric latitude, beta (in degrees)
    let beta_degrees = -b_degrees;

    // 3.3. Calculate the nutation in longitude and obliquity
    let x_terms = calculate_nutation_terms(jce);
    let delta_psi_epsilon = calculate_delta_psi_epsilon(jce, &x_terms);

    // 3.4. Calculate the true obliquity of the ecliptic, epsilon (in degrees)
    let epsilon_degrees =
        calculate_true_obliquity_of_ecliptic(&jd, delta_psi_epsilon.delta_epsilon);

    // 3.5. Calculate the aberration correction, delta_tau (in degrees)
    let delta_tau = ABERRATION_CONSTANT / (SECONDS_PER_HOUR * r);

    // 3.6. Calculate the apparent sun longitude, lambda (in degrees)
    let lambda_degrees = theta_degrees + delta_psi_epsilon.delta_psi + delta_tau;

    // 3.7. Calculate the apparent sidereal time at Greenwich at any given time, nu (in degrees)
    let nu_degrees = calculate_apparent_sidereal_time_at_greenwich(
        &jd,
        delta_psi_epsilon.delta_psi,
        epsilon_degrees,
    );

    // 3.8.1. Calculate the geocentric sun right ascension, alpha (in degrees)
    let beta = degrees_to_radians(beta_degrees);
    let epsilon = degrees_to_radians(epsilon_degrees);
    let lambda = degrees_to_radians(lambda_degrees);
    let alpha_degrees = calculate_geocentric_sun_right_ascension(beta, epsilon, lambda);

    // 3.8.2. Calculate the geocentric sun declination, delta (in degrees)
    let delta_degrees =
        radians_to_degrees(calculate_geocentric_sun_declination(beta, epsilon, lambda));

    Ok(SpaTimeDependent {
        theta_degrees,
        beta_degrees,
        r,
        delta_psi: delta_psi_epsilon.delta_psi,
        epsilon_degrees,
        lambda_degrees,
        nu_degrees,
        alpha_degrees,
        delta_degrees,
    })
}

/// Complete SPA calculation using pre-computed time-dependent parts (steps 12+).
///
/// This function completes the SPA calculation using cached intermediate values
/// from [`spa_time_dependent_parts`]. Used together, these provide significant
/// speedup for coordinate sweeps with unchanged accuracy.
///
/// # Arguments
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `elevation` - Observer elevation above sea level in meters
/// * `refraction` - Optional atmospheric refraction correction
/// * `time_dependent` - Pre-computed time-dependent calculations from [`spa_time_dependent_parts`]
///
/// # Returns
/// Solar position or error
///
/// # Errors
/// Returns error for invalid coordinates (latitude outside ±90°, longitude outside ±180°)
///
/// # Example
/// ```rust
/// use solar_positioning::{spa, RefractionCorrection};
/// use chrono::{DateTime, Utc};
///
/// let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
/// let time_parts = spa::spa_time_dependent_parts(datetime, 69.0).unwrap();
///
/// let position = spa::spa_with_time_dependent_parts(
///     37.7749,   // San Francisco latitude
///     -122.4194, // San Francisco longitude
///     0.0,       // elevation (meters)
///     Some(RefractionCorrection::standard()),
///     &time_parts
/// ).unwrap();
///
/// println!("Azimuth: {:.3}°", position.azimuth());
/// ```
pub fn spa_with_time_dependent_parts(
    latitude: f64,
    longitude: f64,
    elevation: f64,
    refraction: Option<RefractionCorrection>,
    time_dependent: &SpaTimeDependent,
) -> Result<SolarPosition> {
    check_coordinates(latitude, longitude)?;

    // 3.9. Calculate the observer local hour angle, H (in degrees)
    // Use pre-computed apparent sidereal time from time_dependent parts
    let nu_degrees = time_dependent.nu_degrees;

    // Use pre-computed geocentric sun right ascension and declination
    let h_degrees =
        normalize_degrees_0_to_360(nu_degrees + longitude - time_dependent.alpha_degrees);
    let h = degrees_to_radians(h_degrees);

    // 3.10-3.11. Calculate the topocentric sun coordinates
    let xi_degrees = 8.794 / (3600.0 * time_dependent.r);
    let xi = degrees_to_radians(xi_degrees);
    let phi = degrees_to_radians(latitude);
    let delta = degrees_to_radians(time_dependent.delta_degrees);

    let u = atan(EARTH_FLATTENING_FACTOR * tan(phi));
    let y = mul_add(
        EARTH_FLATTENING_FACTOR,
        sin(u),
        (elevation / EARTH_RADIUS_METERS) * sin(phi),
    );
    let x = mul_add(elevation / EARTH_RADIUS_METERS, cos(phi), cos(u));

    let delta_alpha_prime_degrees = radians_to_degrees(atan2(
        -x * sin(xi) * sin(h),
        mul_add(x * sin(xi), -cos(h), cos(delta)),
    ));

    let delta_prime = radians_to_degrees(atan2(
        mul_add(y, -sin(xi), sin(delta)) * cos(degrees_to_radians(delta_alpha_prime_degrees)),
        mul_add(x * sin(xi), -cos(h), cos(delta)),
    ));

    // 3.12. Calculate the topocentric local hour angle, H' (in degrees)
    let h_prime_degrees = h_degrees - delta_alpha_prime_degrees;

    // 3.13. Calculate the topocentric zenith and azimuth angles
    let zenith_angle = radians_to_degrees(acos(mul_add(
        sin(degrees_to_radians(latitude)),
        sin(degrees_to_radians(delta_prime)),
        cos(degrees_to_radians(latitude))
            * cos(degrees_to_radians(delta_prime))
            * cos(degrees_to_radians(h_prime_degrees)),
    )));

    // 3.14. Calculate the topocentric azimuth angle
    let azimuth = normalize_degrees_0_to_360(
        180.0
            + radians_to_degrees(atan2(
                sin(degrees_to_radians(h_prime_degrees)),
                cos(degrees_to_radians(h_prime_degrees)) * sin(degrees_to_radians(latitude))
                    - tan(degrees_to_radians(delta_prime)) * cos(degrees_to_radians(latitude)),
            )),
    );

    // Apply atmospheric refraction if requested
    let elevation_angle = 90.0 - zenith_angle;
    let final_zenith = refraction.map_or(zenith_angle, |correction| {
        if elevation_angle > SUNRISE_SUNSET_ANGLE {
            let pressure = correction.pressure();
            let temperature = correction.temperature();
            // Apply refraction correction following the same pattern as calculate_topocentric_zenith_angle
            zenith_angle
                - (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 1.02
                    / (60.0
                        * tan(degrees_to_radians(
                            elevation_angle + 10.3 / (elevation_angle + 5.11),
                        )))
        } else {
            zenith_angle
        }
    });

    SolarPosition::new(azimuth, final_zenith)
}

#[cfg(all(test, feature = "chrono", feature = "std"))]
mod tests {
    use super::*;
    use chrono::{DateTime, FixedOffset};
    use std::collections::HashSet;

    #[test]
    fn test_spa_basic_functionality() {
        let datetime = "2023-06-21T12:00:00Z"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();

        let result = solar_position(
            datetime,
            37.7749, // San Francisco
            -122.4194,
            0.0,
            69.0,
            Some(RefractionCorrection::new(1013.25, 15.0).unwrap()),
        );

        assert!(result.is_ok());
        let position = result.unwrap();
        assert!(position.azimuth() >= 0.0 && position.azimuth() <= 360.0);
        assert!(position.zenith_angle() >= 0.0 && position.zenith_angle() <= 180.0);
    }

    #[test]
    fn test_sunrise_sunset_multiple() {
        let datetime = "2023-06-21T12:00:00Z"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();
        let horizons = [
            Horizon::SunriseSunset,
            Horizon::CivilTwilight,
            Horizon::NauticalTwilight,
        ];

        let results: Result<Vec<_>> = sunrise_sunset_multiple(
            datetime,
            37.7749,   // San Francisco latitude
            -122.4194, // San Francisco longitude
            69.0,      // deltaT (seconds)
            horizons.iter().copied(),
        )
        .collect();

        let results = results.unwrap();

        // Should have results for all requested horizons
        assert_eq!(results.len(), 3);

        // Check that we have all expected horizons
        let returned_horizons: HashSet<_> = results.iter().map(|(h, _)| *h).collect();
        for expected_horizon in horizons {
            assert!(returned_horizons.contains(&expected_horizon));
        }

        // Compare with individual calls to ensure consistency
        for (horizon, bulk_result) in &results {
            let individual_result =
                sunrise_sunset_for_horizon(datetime, 37.7749, -122.4194, 69.0, *horizon).unwrap();

            // Results should be identical
            match (&individual_result, bulk_result) {
                (
                    crate::SunriseResult::RegularDay {
                        sunrise: s1,
                        transit: t1,
                        sunset: ss1,
                    },
                    crate::SunriseResult::RegularDay {
                        sunrise: s2,
                        transit: t2,
                        sunset: ss2,
                    },
                ) => {
                    assert_eq!(s1, s2);
                    assert_eq!(t1, t2);
                    assert_eq!(ss1, ss2);
                }
                (
                    crate::SunriseResult::AllDay { transit: t1 },
                    crate::SunriseResult::AllDay { transit: t2 },
                )
                | (
                    crate::SunriseResult::AllNight { transit: t1 },
                    crate::SunriseResult::AllNight { transit: t2 },
                ) => {
                    assert_eq!(t1, t2);
                }
                _ => panic!("Bulk and individual results differ in type for {horizon:?}"),
            }
        }
    }

    #[test]
    fn test_spa_no_refraction() {
        let datetime = "2023-06-21T12:00:00Z"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();

        let result = solar_position(datetime, 37.7749, -122.4194, 0.0, 69.0, None);

        assert!(result.is_ok());
        let position = result.unwrap();
        assert!(position.azimuth() >= 0.0 && position.azimuth() <= 360.0);
        assert!(position.zenith_angle() >= 0.0 && position.zenith_angle() <= 180.0);
    }

    #[test]
    fn test_spa_coordinate_validation() {
        let datetime = "2023-06-21T12:00:00Z"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();

        // Invalid latitude
        assert!(solar_position(
            datetime,
            95.0,
            0.0,
            0.0,
            0.0,
            Some(RefractionCorrection::new(1013.25, 15.0).unwrap())
        )
        .is_err());

        // Invalid longitude
        assert!(solar_position(
            datetime,
            0.0,
            185.0,
            0.0,
            0.0,
            Some(RefractionCorrection::new(1013.25, 15.0).unwrap())
        )
        .is_err());
    }

    #[test]
    fn test_sunrise_sunset_basic() {
        let date = "2023-06-21T00:00:00Z"
            .parse::<DateTime<FixedOffset>>()
            .unwrap();

        let result = sunrise_sunset(date, 37.7749, -122.4194, 69.0, -0.833);
        assert!(result.is_ok());

        let result =
            sunrise_sunset_for_horizon(date, 37.7749, -122.4194, 69.0, Horizon::SunriseSunset);
        assert!(result.is_ok());
    }

    #[test]
    fn test_horizon_enum() {
        assert_eq!(Horizon::SunriseSunset.elevation_angle(), -0.83337);
        assert_eq!(Horizon::CivilTwilight.elevation_angle(), -6.0);
        assert_eq!(Horizon::Custom(-10.5).elevation_angle(), -10.5);
    }
}
