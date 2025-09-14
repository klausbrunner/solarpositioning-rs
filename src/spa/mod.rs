//! SPA (Solar Position Algorithm) implementation.
//!
//! High-accuracy solar positioning based on the NREL algorithm by Reda & Andreas (2003).
//! Provides uncertainties of ±0.0003 degrees for years -2000 to 6000.

#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::unreadable_literal)]

use crate::error::{check_coordinates, check_refraction_params_usable};
use crate::math::{
    acos, asin, atan, atan2, cos, degrees_to_radians, floor, normalize_degrees_0_to_360,
    polynomial, radians_to_degrees, sin, tan,
};
use crate::time::JulianDate;
use crate::{Error, Horizon, Result, SolarPosition};

pub mod coefficients;
use coefficients::{
    NUTATION_COEFFS, OBLIQUITY_COEFFS, TERMS_B, TERMS_L, TERMS_PE, TERMS_R, TERMS_Y,
};

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
/// * `datetime` - UTC date and time
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `elevation` - Observer elevation in meters above sea level
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `pressure` - Atmospheric pressure in millibars (for refraction correction)
/// * `temperature` - Temperature in °C (for refraction correction)
///
/// # Returns
/// Solar position or error
///
/// # Errors
/// Returns error for invalid coordinates (latitude outside ±90°, longitude outside ±180°)
///
/// # Example
/// ```rust
/// use solar_positioning::spa;
/// use chrono::{DateTime, Utc};
///
/// let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
/// let position = spa::solar_position(
///     datetime,
///     37.7749,     // San Francisco latitude
///     -122.4194,   // San Francisco longitude
///     0.0,         // elevation (meters)
///     69.0,        // deltaT (seconds)
///     1013.25,     // pressure (millibars)
///     15.0         // temperature (°C)
/// ).unwrap();
///
/// println!("Azimuth: {:.3}°", position.azimuth());
/// println!("Elevation: {:.3}°", position.elevation_angle());
/// ```
pub fn solar_position(
    datetime: chrono::DateTime<chrono::Utc>,
    latitude: f64,
    longitude: f64,
    elevation: f64,
    delta_t: f64,
    pressure: f64,
    temperature: f64,
) -> Result<SolarPosition> {
    check_coordinates(latitude, longitude)?;

    // Convert datetime to Julian date components
    let jd = JulianDate::from_datetime(datetime, delta_t)?;

    let jme = jd.julian_ephemeris_millennium();
    let jce = jd.julian_ephemeris_century();

    // Step 1: Calculate Earth heliocentric longitude, L
    let l_terms = calculate_lbr_terms(jme, TERMS_L);
    let l_degrees =
        normalize_degrees_0_to_360(radians_to_degrees(calculate_lbr_polynomial(jme, &l_terms)));

    // Step 2: Calculate Earth heliocentric latitude, B
    let b_terms = calculate_lbr_terms(jme, TERMS_B);
    let b_degrees =
        normalize_degrees_0_to_360(radians_to_degrees(calculate_lbr_polynomial(jme, &b_terms)));

    // Step 3: Calculate Earth radius vector, R
    let r_terms = calculate_lbr_terms(jme, TERMS_R);
    let r = calculate_lbr_polynomial(jme, &r_terms);

    if r == 0.0 {
        return Err(Error::computation_error("Earth radius vector is zero"));
    }

    // Step 4: Calculate geocentric longitude and latitude
    let theta_degrees = normalize_degrees_0_to_360(l_degrees + 180.0);
    let beta_degrees = -b_degrees;
    let beta = degrees_to_radians(beta_degrees);

    // Step 5: Calculate nutation
    let x_terms = calculate_nutation_terms(jce);
    let delta_psi_epsilon = calculate_delta_psi_epsilon(jce, &x_terms);

    // Step 6: Calculate true obliquity of the ecliptic
    let epsilon_degrees =
        calculate_true_obliquity_of_ecliptic(&jd, delta_psi_epsilon.delta_epsilon);
    let epsilon = degrees_to_radians(epsilon_degrees);

    // Step 7: Calculate aberration correction
    let delta_tau = ABERRATION_CONSTANT / (SECONDS_PER_HOUR * r);

    // Step 8: Calculate apparent sun longitude
    let lambda_degrees = theta_degrees + delta_psi_epsilon.delta_psi + delta_tau;
    let lambda = degrees_to_radians(lambda_degrees);

    // Step 9: Calculate apparent sidereal time at Greenwich
    let nu_degrees = calculate_apparent_sidereal_time_at_greenwich(
        &jd,
        delta_psi_epsilon.delta_psi,
        epsilon_degrees,
    );

    // Step 10: Calculate geocentric sun right ascension
    let alpha_degrees = calculate_geocentric_sun_right_ascension(beta, epsilon, lambda);

    // Step 11: Calculate geocentric sun declination
    let delta_degrees =
        radians_to_degrees(calculate_geocentric_sun_declination(beta, epsilon, lambda));

    // Step 12: Calculate observer local hour angle
    let h_degrees = normalize_degrees_0_to_360(nu_degrees + longitude - alpha_degrees);
    let h = degrees_to_radians(h_degrees);

    // Step 13: Calculate topocentric sun right ascension and declination
    let xi_degrees = 8.794 / (3600.0 * r);
    let xi = degrees_to_radians(xi_degrees);
    let phi = degrees_to_radians(latitude);
    let delta = degrees_to_radians(delta_degrees);

    let u = atan(EARTH_FLATTENING_FACTOR * tan(phi));
    let x = cos(u) + elevation * cos(phi) / EARTH_RADIUS_METERS;
    let y = EARTH_FLATTENING_FACTOR.mul_add(sin(u), (elevation * sin(phi)) / EARTH_RADIUS_METERS);

    let x1 = (x * sin(xi)).mul_add(-cos(h), cos(delta));
    let delta_alpha_degrees = radians_to_degrees(atan2(-x * sin(xi) * sin(h), x1));
    let delta_prime = atan2(
        y.mul_add(-sin(xi), sin(delta)) * cos(degrees_to_radians(delta_alpha_degrees)),
        x1,
    );

    // Step 14: Calculate topocentric local hour angle
    let h_prime_degrees = h_degrees - delta_alpha_degrees;
    let h_prime = degrees_to_radians(h_prime_degrees);

    // Step 15: Calculate topocentric solar position
    calculate_topocentric_solar_position(pressure, temperature, phi, delta_prime, h_prime)
}

/// Calculate solar position using the SPA algorithm without refraction correction.
///
/// This is identical to the main function but doesn't apply atmospheric refraction correction.
///
/// # Errors
/// Returns error for invalid coordinates (latitude outside ±90°, longitude outside ±180°)
pub fn solar_position_no_refraction(
    datetime: chrono::DateTime<chrono::Utc>,
    latitude: f64,
    longitude: f64,
    elevation: f64,
    delta_t: f64,
) -> Result<SolarPosition> {
    solar_position(
        datetime,
        latitude,
        longitude,
        elevation,
        delta_t,
        f64::NAN,
        f64::NAN,
    )
}

/// Represents nutation corrections for longitude and obliquity.
#[derive(Debug, Clone, Copy)]
struct DeltaPsiEpsilon {
    delta_psi: f64,
    delta_epsilon: f64,
}

/// Calculate L, B, R terms from the coefficient tables.
fn calculate_lbr_terms(jme: f64, term_coeffs: &[&[&[f64; 3]]]) -> Vec<f64> {
    let mut lbr_terms = vec![0.0; term_coeffs.len()];

    for (i, term_set) in term_coeffs.iter().enumerate() {
        let mut lbr_sum = 0.0;
        for term in *term_set {
            lbr_sum += term[0] * cos(term[2].mul_add(jme, term[1]));
        }
        lbr_terms[i] = lbr_sum;
    }

    lbr_terms
}

/// Calculate L, B, or R polynomial from the terms.
fn calculate_lbr_polynomial(jme: f64, terms: &[f64]) -> f64 {
    polynomial(terms, jme) / 1e8
}

/// Calculate nutation terms (X values).
fn calculate_nutation_terms(jce: f64) -> Vec<f64> {
    NUTATION_COEFFS
        .iter()
        .map(|coeffs| polynomial(coeffs, jce))
        .collect()
}

/// Calculate nutation in longitude and obliquity.
fn calculate_delta_psi_epsilon(jce: f64, x: &[f64]) -> DeltaPsiEpsilon {
    let mut delta_psi = 0.0;
    let mut delta_epsilon = 0.0;

    for (i, pe_term) in TERMS_PE.iter().enumerate() {
        let xj_yterm_sum = degrees_to_radians(calculate_xj_yterm_sum(i, x));

        // Use Math.fma equivalent: a * b + c
        let delta_psi_contrib = pe_term[1].mul_add(jce, pe_term[0]) * sin(xj_yterm_sum);
        let delta_epsilon_contrib = pe_term[3].mul_add(jce, pe_term[2]) * cos(xj_yterm_sum);

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
    let nu0_degrees = normalize_degrees_0_to_360(jd.julian_century().powi(2).mul_add(
        0.000387933 - jd.julian_century() / 38710000.0,
        360.98564736629f64.mul_add(jd.julian_date() - 2451545.0, 280.46061837),
    ));

    delta_psi.mul_add(cos(degrees_to_radians(epsilon_degrees)), nu0_degrees)
}

/// Calculate geocentric sun right ascension.
fn calculate_geocentric_sun_right_ascension(
    beta_rad: f64,
    epsilon_rad: f64,
    lambda_rad: f64,
) -> f64 {
    let alpha = atan2(
        sin(lambda_rad).mul_add(cos(epsilon_rad), -(tan(beta_rad) * sin(epsilon_rad))),
        cos(lambda_rad),
    );
    normalize_degrees_0_to_360(radians_to_degrees(alpha))
}

/// Calculate geocentric sun declination.
fn calculate_geocentric_sun_declination(beta_rad: f64, epsilon_rad: f64, lambda_rad: f64) -> f64 {
    asin(sin(beta_rad).mul_add(
        cos(epsilon_rad),
        cos(beta_rad) * sin(epsilon_rad) * sin(lambda_rad),
    ))
}

/// Calculate topocentric solar position with optional refraction correction.
fn calculate_topocentric_solar_position(
    pressure: f64,
    temperature: f64,
    phi: f64,
    delta_prime: f64,
    h_prime: f64,
) -> Result<SolarPosition> {
    // Calculate topocentric zenith angle
    let sin_phi = sin(phi);
    let cos_phi = cos(phi);
    let cos_h_prime = cos(h_prime);

    let e_zero = asin(sin_phi.mul_add(sin(delta_prime), cos_phi * cos(delta_prime) * cos_h_prime));
    let topocentric_zenith_angle =
        calculate_topocentric_zenith_angle(pressure, temperature, e_zero);

    // Calculate topocentric azimuth angle
    let gamma = atan2(
        sin(h_prime),
        cos_h_prime.mul_add(sin_phi, -(tan(delta_prime) * cos_phi)),
    );
    let gamma_degrees = normalize_degrees_0_to_360(radians_to_degrees(gamma));
    let topocentric_azimuth_angle = normalize_degrees_0_to_360(gamma_degrees + 180.0);

    SolarPosition::new(topocentric_azimuth_angle, topocentric_zenith_angle)
}

/// Calculate topocentric zenith angle with optional refraction correction.
fn calculate_topocentric_zenith_angle(pressure: f64, temperature: f64, e_zero: f64) -> f64 {
    let e_zero_degrees = radians_to_degrees(e_zero);

    // Apply refraction correction if parameters are usable and sun is visible
    let do_correct = check_refraction_params_usable(pressure, temperature)
        && e_zero_degrees > SUNRISE_SUNSET_ANGLE;

    if do_correct {
        90.0 - e_zero_degrees
            - (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 1.02
                / (60.0
                    * tan(degrees_to_radians(
                        e_zero_degrees + 10.3 / (e_zero_degrees + 5.11),
                    )))
    } else {
        90.0 - e_zero_degrees
    }
}

/// Calculate sunrise, solar transit, and sunset times using the SPA algorithm.
///
/// This follows the NREL SPA algorithm (Reda & Andreas 2003) for calculating
/// sunrise, transit (solar noon), and sunset times with high accuracy.
///
/// # Arguments
/// * `date` - UTC date for calculations (time is ignored, calculations are for the entire day)
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
/// May panic if date conversion fails (should not occur for valid UTC dates)
///
/// # Example
/// ```rust
/// use solar_positioning::spa;
/// use chrono::{DateTime, Utc, NaiveDate};
///
/// let date = NaiveDate::from_ymd_opt(2023, 6, 21).unwrap()
///     .and_hms_opt(0, 0, 0).unwrap()
///     .and_utc();
/// let result = spa::sunrise_sunset(
///     date,
///     37.7749,   // San Francisco latitude
///     -122.4194, // San Francisco longitude
///     69.0,      // deltaT (seconds)
///     -0.833     // standard sunrise/sunset angle
/// ).unwrap();
/// ```
pub fn sunrise_sunset(
    date: chrono::DateTime<chrono::Utc>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    elevation_angle: f64,
) -> Result<crate::SunriseResult> {
    check_coordinates(latitude, longitude)?;

    // Create Julian date for the day
    let day_start = date.date_naive().and_hms_opt(0, 0, 0).unwrap().and_utc();

    // Implement the complete Java algorithm in one place for accuracy
    calculate_sunrise_sunset_spa_algorithm(day_start, latitude, longitude, delta_t, elevation_angle)
}

/// Calculate sunrise/sunset times using NREL SPA algorithm Appendix A.2
fn calculate_sunrise_sunset_spa_algorithm(
    day: chrono::DateTime<chrono::Utc>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    elevation_angle: f64,
) -> Result<crate::SunriseResult> {
    let day_start = day.date_naive().and_hms_opt(0, 0, 0).unwrap().and_utc();

    // A.2.1. Calculate the apparent sidereal time at Greenwich at 0 UT
    let (nu_degrees, delta_psi_epsilon, epsilon_degrees) =
        calculate_sidereal_time_and_nutation(day_start);

    // A.2.2. Calculate alpha/delta for day before, same day, next day
    let alpha_deltas =
        calculate_alpha_deltas_for_three_days(day_start, delta_psi_epsilon, epsilon_degrees)?;

    // Calculate initial transit time and check for polar conditions
    let m0 = (alpha_deltas[1].alpha - longitude - nu_degrees) / 360.0;
    if let Some(polar_result) =
        check_polar_conditions(day, m0, latitude, elevation_angle, alpha_deltas[1].delta)
    {
        return Ok(polar_result);
    }

    // Calculate approximate times and apply corrections
    let (m_values, _h0_degrees) =
        calculate_approximate_times(m0, latitude, elevation_angle, alpha_deltas[1].delta);

    // Apply final corrections to get accurate times
    Ok(calculate_final_times(
        day,
        m_values,
        nu_degrees,
        delta_t,
        latitude,
        longitude,
        elevation_angle,
        alpha_deltas,
    ))
}

/// A.2.1. Calculate apparent sidereal time and nutation parameters
fn calculate_sidereal_time_and_nutation(
    day_start: chrono::DateTime<chrono::Utc>,
) -> (f64, DeltaPsiEpsilon, f64) {
    let jd = JulianDate::from_datetime(day_start, 0.0).unwrap();
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

/// A.2.2. Calculate alpha/delta for day before, same day, and next day
fn calculate_alpha_deltas_for_three_days(
    day_start: chrono::DateTime<chrono::Utc>,
    delta_psi_epsilon: DeltaPsiEpsilon,
    epsilon_degrees: f64,
) -> Result<[AlphaDelta; 3]> {
    let mut alpha_deltas = [AlphaDelta {
        alpha: 0.0,
        delta: 0.0,
    }; 3];
    for (i, alpha_delta) in alpha_deltas.iter_mut().enumerate() {
        let current_jd = JulianDate::from_datetime(day_start, 0.0)?.add_days((i as f64) - 1.0);
        let current_jme = current_jd.julian_ephemeris_millennium();
        let ad = calculate_alpha_delta(current_jme, delta_psi_epsilon.delta_psi, epsilon_degrees);
        *alpha_delta = ad;
    }
    Ok(alpha_deltas)
}

/// Check for polar day/night conditions and return early if they apply
fn check_polar_conditions(
    day: chrono::DateTime<chrono::Utc>,
    m0: f64,
    latitude: f64,
    elevation_angle: f64,
    delta1: f64,
) -> Option<crate::SunriseResult> {
    let phi = degrees_to_radians(latitude);
    let delta1_rad = degrees_to_radians(delta1);
    let elevation_rad = degrees_to_radians(elevation_angle);

    let acos_arg =
        sin(phi).mul_add(-sin(delta1_rad), sin(elevation_rad)) / (cos(phi) * cos(delta1_rad));

    if acos_arg < -1.0 {
        let transit = add_fraction_of_day(day, m0);
        Some(crate::SunriseResult::AllDay { transit })
    } else if acos_arg > 1.0 {
        let transit = add_fraction_of_day(day, m0);
        Some(crate::SunriseResult::AllNight { transit })
    } else {
        None
    }
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
        sin(phi).mul_add(-sin(delta1_rad), sin(elevation_rad)) / (cos(phi) * cos(delta1_rad));
    let h0 = acos(acos_arg);
    let h0_degrees = radians_to_degrees(h0).min(180.0);

    let mut m = [0.0; 3];
    m[0] = normalize_to_unit_range(m0);
    m[1] = normalize_to_unit_range(m0 - h0_degrees / 360.0);
    m[2] = normalize_to_unit_range(m0 + h0_degrees / 360.0);

    (m, h0_degrees)
}

/// A.2.8-15. Calculate final accurate times using corrections
fn calculate_final_times(
    day: chrono::DateTime<chrono::Utc>,
    m_values: [f64; 3],
    nu_degrees: f64,
    delta_t: f64,
    latitude: f64,
    longitude: f64,
    elevation_angle: f64,
    alpha_deltas: [AlphaDelta; 3],
) -> crate::SunriseResult {
    // A.2.8. Calculate sidereal times
    let mut nu = [0.0; 3];
    for i in 0..3 {
        nu[i] = 360.985647f64.mul_add(m_values[i], nu_degrees);
    }

    // A.2.9. Calculate terms with deltaT correction
    let mut n = [0.0; 3];
    for i in 0..3 {
        n[i] = m_values[i] + delta_t / 86400.0;
    }

    // A.2.10. Calculate α'i and δ'i using quadratic interpolation
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
        h[i] = radians_to_degrees(asin(sin(phi).mul_add(
            sin(delta_prime_rad),
            cos(phi) * cos(delta_prime_rad) * cos(degrees_to_radians(h_prime[i])),
        )));
    }

    // A.2.13-15. Calculate final times
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

    let sunrise = add_fraction_of_day(day, r);
    let transit = add_fraction_of_day(day, t);
    let sunset = add_fraction_of_day(day, s);

    crate::SunriseResult::RegularDay {
        sunrise,
        transit,
        sunset,
    }
}

/// A.2.10. Calculate interpolated alpha/delta values using quadratic interpolation
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
            alpha_deltas[1].alpha + (n[i] * (c.mul_add(n[i], a + b))) / 2.0;
        alpha_delta_primes[i].delta =
            alpha_deltas[1].delta + (n[i] * (c_prime.mul_add(n[i], a_prime + b_prime))) / 2.0;
    }
    alpha_delta_primes
}

#[derive(Debug, Clone, Copy)]
struct AlphaDelta {
    alpha: f64,
    delta: f64,
}

/// Calculate sunrise, solar transit, and sunset times for a specific horizon type.
///
/// This is a convenience function that uses predefined elevation angles for common
/// sunrise/twilight calculations.
///
/// # Arguments
/// * `date` - UTC date for calculations
/// * `latitude` - Observer latitude in degrees
/// * `longitude` - Observer longitude in degrees
/// * `delta_t` - ΔT in seconds
/// * `horizon` - Horizon type (sunrise/sunset, civil twilight, etc.)
///
/// # Errors
/// Returns error for invalid coordinates, dates, or computation failures.
///
/// # Example
/// ```rust
/// use solar_positioning::{spa, Horizon};
/// use chrono::{NaiveDate};
///
/// let date = NaiveDate::from_ymd_opt(2023, 6, 21).unwrap()
///     .and_hms_opt(0, 0, 0).unwrap()
///     .and_utc();
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
pub fn sunrise_sunset_for_horizon(
    date: chrono::DateTime<chrono::Utc>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    horizon: Horizon,
) -> Result<crate::SunriseResult> {
    sunrise_sunset(
        date,
        latitude,
        longitude,
        delta_t,
        horizon.elevation_angle(),
    )
}

/// Calculate alpha (right ascension) and delta (declination) for a given JME using full SPA algorithm
fn calculate_alpha_delta(jme: f64, delta_psi: f64, epsilon_degrees: f64) -> AlphaDelta {
    // Follow Java calculateAlphaDelta exactly

    // calculate Earth heliocentric latitude, B
    let b_terms = calculate_lbr_terms(jme, TERMS_B);
    let b_degrees =
        normalize_degrees_0_to_360(radians_to_degrees(calculate_lbr_polynomial(jme, &b_terms)));

    // calculate Earth radius vector, R
    let r_terms = calculate_lbr_terms(jme, TERMS_R);
    let r = calculate_lbr_polynomial(jme, &r_terms);
    assert!(r != 0.0);

    // calculate Earth heliocentric longitude, L
    let l_terms = calculate_lbr_terms(jme, TERMS_L);
    let l_degrees =
        normalize_degrees_0_to_360(radians_to_degrees(calculate_lbr_polynomial(jme, &l_terms)));

    // calculate geocentric longitude, theta
    let theta_degrees = normalize_degrees_0_to_360(l_degrees + 180.0);

    // calculate geocentric latitude, beta
    let beta_degrees = -b_degrees;
    let beta = degrees_to_radians(beta_degrees);
    let epsilon = degrees_to_radians(epsilon_degrees);

    // calculate aberration correction
    let delta_tau = ABERRATION_CONSTANT / (SECONDS_PER_HOUR * r);

    // calculate the apparent sun longitude
    let lambda_degrees = theta_degrees + delta_psi + delta_tau;
    let lambda = degrees_to_radians(lambda_degrees);

    // Calculate the geocentric sun right ascension and declination (like Java)
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
    let limited = divided - divided.floor();
    if limited < 0.0 {
        limited + 1.0
    } else {
        limited
    }
}

/// Add a fraction of a day to a date
fn add_fraction_of_day(
    day: chrono::DateTime<chrono::Utc>,
    fraction: f64,
) -> chrono::DateTime<chrono::Utc> {
    let millis = (fraction * 24.0 * 60.0 * 60.0 * 1000.0) as i64;
    day + chrono::Duration::milliseconds(millis)
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
/// * `date` - UTC date and time
/// * `latitude` - Observer latitude in degrees (-90 to +90)
/// * `longitude` - Observer longitude in degrees (-180 to +180)
/// * `delta_t` - ΔT in seconds (difference between TT and UT1)
/// * `horizons` - Iterator of horizon types to calculate
///
/// # Returns
/// Iterator over `Result<(Horizon, SunriseResult)>`
///
/// # Panics
/// May panic if date cannot be converted to start of day (unlikely with valid dates).
///
/// # Example
/// ```rust
/// use solar_positioning::{spa, Horizon};
/// use chrono::{DateTime, Utc};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
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
pub fn sunrise_sunset_multiple<H>(
    date: chrono::DateTime<chrono::Utc>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    horizons: H,
) -> impl Iterator<Item = Result<(Horizon, crate::SunriseResult)>>
where
    H: IntoIterator<Item = Horizon>,
{
    // Pre-calculate common values once for efficiency
    let precomputed = (|| -> Result<_> {
        check_coordinates(latitude, longitude)?;

        let day_start = date.date_naive().and_hms_opt(0, 0, 0).unwrap().and_utc();
        let (nu_degrees, delta_psi_epsilon, epsilon_degrees) =
            calculate_sidereal_time_and_nutation(day_start);
        let alpha_deltas =
            calculate_alpha_deltas_for_three_days(day_start, delta_psi_epsilon, epsilon_degrees)?;

        Ok((nu_degrees, alpha_deltas))
    })();

    horizons.into_iter().map(move |horizon| match &precomputed {
        Ok((nu_degrees, alpha_deltas)) => {
            let result = calculate_sunrise_sunset_spa_algorithm_with_precomputed(
                date,
                latitude,
                longitude,
                delta_t,
                horizon.elevation_angle(),
                *nu_degrees,
                *alpha_deltas,
            )?;
            Ok((horizon, result))
        }
        Err(e) => Err(e.clone()),
    })
}

/// Internal helper that uses precomputed values for efficiency
fn calculate_sunrise_sunset_spa_algorithm_with_precomputed(
    day: chrono::DateTime<chrono::Utc>,
    latitude: f64,
    longitude: f64,
    delta_t: f64,
    elevation_angle: f64,
    nu_degrees: f64,
    alpha_deltas: [AlphaDelta; 3],
) -> Result<crate::SunriseResult> {
    // Use the same day_start as the individual calculation to ensure consistency
    let day_start = day.date_naive().and_hms_opt(0, 0, 0).unwrap().and_utc();

    // Calculate initial transit time and check for polar conditions
    let m0 = (alpha_deltas[1].alpha - longitude - nu_degrees) / 360.0;
    if let Some(polar_result) = check_polar_conditions(
        day_start,
        m0,
        latitude,
        elevation_angle,
        alpha_deltas[1].delta,
    ) {
        return Ok(polar_result);
    }

    // Calculate approximate times and apply corrections
    let (m_values, _h0_degrees) =
        calculate_approximate_times(m0, latitude, elevation_angle, alpha_deltas[1].delta);

    // Apply final corrections to get accurate times
    Ok(calculate_final_times(
        day_start,
        m_values,
        nu_degrees,
        delta_t,
        latitude,
        longitude,
        elevation_angle,
        alpha_deltas,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{DateTime, Utc};

    #[test]
    fn test_spa_basic_functionality() {
        let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();

        let result = solar_position(
            datetime, 37.7749, // San Francisco
            -122.4194, 0.0, 69.0, 1013.25, 15.0,
        );

        assert!(result.is_ok());
        let position = result.unwrap();
        assert!(position.azimuth() >= 0.0 && position.azimuth() <= 360.0);
        assert!(position.zenith_angle() >= 0.0 && position.zenith_angle() <= 180.0);
    }

    #[test]
    fn test_sunrise_sunset_multiple() {
        let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();
        let horizons = [
            Horizon::SunriseSunset,
            Horizon::CivilTwilight,
            Horizon::NauticalTwilight,
        ];

        let results: std::result::Result<Vec<_>, _> = sunrise_sunset_multiple(
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
        let returned_horizons: std::collections::HashSet<_> =
            results.iter().map(|(h, _)| *h).collect();
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
        let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();

        let result = solar_position_no_refraction(datetime, 37.7749, -122.4194, 0.0, 69.0);

        assert!(result.is_ok());
        let position = result.unwrap();
        assert!(position.azimuth() >= 0.0 && position.azimuth() <= 360.0);
        assert!(position.zenith_angle() >= 0.0 && position.zenith_angle() <= 180.0);
    }

    #[test]
    fn test_spa_coordinate_validation() {
        let datetime = "2023-06-21T12:00:00Z".parse::<DateTime<Utc>>().unwrap();

        // Invalid latitude
        assert!(solar_position(datetime, 95.0, 0.0, 0.0, 0.0, 1013.25, 15.0).is_err());

        // Invalid longitude
        assert!(solar_position(datetime, 0.0, 185.0, 0.0, 0.0, 1013.25, 15.0).is_err());
    }

    #[test]
    fn test_sunrise_sunset_basic() {
        let date = "2023-06-21T00:00:00Z".parse::<DateTime<Utc>>().unwrap();

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
