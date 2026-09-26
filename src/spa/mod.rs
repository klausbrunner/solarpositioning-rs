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
    acos, asin, atan, atan2, cos, degrees_to_radians, mul_add, normalize_degrees_0_to_360,
    polynomial, powi, radians_to_degrees, sin, sin_cos, sqrt, tan,
};
use crate::time::JulianDate;
use crate::{Location, RefractionCorrection, Result, SolarPosition, events::EventPosition};

mod coefficients;
use coefficients::{
    NUTATION_COEFFS, OBLIQUITY_COEFFS, TERMS_B, TERMS_L, TERMS_PE, TERMS_R, TERMS_Y,
};

/// Aberration constant in arcseconds.
const ABERRATION_CONSTANT: f64 = -20.4898;

/// Earth flattening factor (WGS84).
const EARTH_FLATTENING_FACTOR: f64 = 0.99664719;

/// Earth radius in meters (WGS84).
const EARTH_RADIUS_METERS: f64 = 6378140.0;

/// Seconds per hour conversion factor.
const SECONDS_PER_HOUR: f64 = 3600.0;

/// Time-dependent SPA values (steps 1-11), independent of observer location.
#[derive(Debug, Clone, Copy)]
pub struct TimeDependent {
    r: f64,
    nu_degrees: f64,
    alpha_degrees: f64,
    delta_degrees: f64,
}

#[derive(Debug, Clone, Copy)]
struct DeltaPsiEpsilon {
    delta_psi: f64,
    delta_epsilon: f64,
}

/// Calculate L, B, or R polynomial from the terms.
fn calculate_lbr_polynomial(jme: f64, term_coeffs: &[&[&[f64; 3]]]) -> f64 {
    let mut term_sums = [0.0; 6];

    for (i, term_set) in term_coeffs.iter().enumerate() {
        let mut sum = 0.0;
        for term in *term_set {
            sum = mul_add(term[0], cos(mul_add(term[2], jme, term[1])), sum);
        }
        term_sums[i] = sum;
    }

    polynomial(&term_sums[..term_coeffs.len()], jme) / 1e8
}

/// Calculate normalized degrees from LBR polynomial
fn lbr_to_normalized_degrees(jme: f64, term_coeffs: &[&[&[f64; 3]]]) -> f64 {
    normalize_degrees_0_to_360(radians_to_degrees(calculate_lbr_polynomial(
        jme,
        term_coeffs,
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
fn calculate_delta_psi_epsilon(jce: f64, x: &[f64; 5]) -> DeltaPsiEpsilon {
    let mut delta_psi = 0.0;
    let mut delta_epsilon = 0.0;

    for (i, pe_term) in TERMS_PE.iter().enumerate() {
        let mut xj_yterm_sum = 0.0;
        for (j, &x_val) in x.iter().enumerate() {
            xj_yterm_sum = mul_add(x_val, f64::from(TERMS_Y[i][j]), xj_yterm_sum);
        }
        let xj_yterm_sum = degrees_to_radians(xj_yterm_sum);

        // Use Math.fma equivalent: a * b + c
        let (sin_sum, cos_sum) = sin_cos(xj_yterm_sum);
        let delta_psi_contrib = mul_add(pe_term[1], jce, pe_term[0]) * sin_sum;
        let delta_epsilon_contrib = mul_add(pe_term[3], jce, pe_term[2]) * cos_sum;

        delta_psi += delta_psi_contrib;
        delta_epsilon += delta_epsilon_contrib;
    }

    DeltaPsiEpsilon {
        delta_psi: delta_psi / 36_000_000.0,
        delta_epsilon: delta_epsilon / 36_000_000.0,
    }
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

/// Calculate geocentric sun right ascension and declination.
fn calculate_geocentric_sun_coordinates(
    beta_rad: f64,
    epsilon_rad: f64,
    lambda_rad: f64,
) -> (f64, f64) {
    let (sin_lambda, cos_lambda) = sin_cos(lambda_rad);
    let (sin_epsilon, cos_epsilon) = sin_cos(epsilon_rad);
    let (sin_beta, cos_beta) = sin_cos(beta_rad);

    let alpha = atan2(
        mul_add(
            sin_lambda,
            cos_epsilon,
            -(sin_beta / cos_beta) * sin_epsilon,
        ),
        cos_lambda,
    );
    let delta = asin(mul_add(
        sin_beta,
        cos_epsilon,
        cos_beta * sin_epsilon * sin_lambda,
    ));
    (
        normalize_degrees_0_to_360(radians_to_degrees(alpha)),
        radians_to_degrees(delta),
    )
}

/// Calculates the astronomical quantities independent of observer location.
pub fn time_dependent(jd: JulianDate) -> TimeDependent {
    let jme = jd.julian_ephemeris_millennium();
    let jce = jd.julian_ephemeris_century();

    // 3.2.2. Calculate the Earth heliocentric longitude, L (in degrees)
    let l_degrees = lbr_to_normalized_degrees(jme, TERMS_L);

    // 3.2.3. Calculate the Earth heliocentric latitude, B (in degrees)
    let b_degrees = lbr_to_normalized_degrees(jme, TERMS_B);

    // 3.2.4. Calculate the Earth radius vector, R (in Astronomical Units, AU)
    let r = calculate_lbr_polynomial(jme, TERMS_R);

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
    let (alpha_degrees, delta_degrees) =
        calculate_geocentric_sun_coordinates(beta, epsilon, lambda);

    TimeDependent {
        r,
        nu_degrees,
        alpha_degrees,
        delta_degrees,
    }
}

/// Completes SPA (steps 12+) for one location using cached time-dependent values.
pub fn solar_position(
    latitude: f64,
    longitude: f64,
    elevation: f64,
    refraction: Option<RefractionCorrection>,
    time_dependent: &TimeDependent,
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
    let sin_xi = sin(xi);
    let (sin_phi, cos_phi) = sin_cos(phi);
    let (sin_delta, cos_delta) = sin_cos(delta);
    let (sin_h, cos_h) = sin_cos(h);

    let u = atan(EARTH_FLATTENING_FACTOR * tan(phi));
    let (sin_u, cos_u) = sin_cos(u);
    let y = mul_add(
        EARTH_FLATTENING_FACTOR,
        sin_u,
        (elevation / EARTH_RADIUS_METERS) * sin_phi,
    );
    let x = mul_add(elevation / EARTH_RADIUS_METERS, cos_phi, cos_u);

    let delta_alpha_prime_degrees = radians_to_degrees(atan2(
        -x * sin_xi * sin_h,
        mul_add(x * sin_xi, -cos_h, cos_delta),
    ));

    let delta_prime_degrees = radians_to_degrees(atan2(
        mul_add(y, -sin_xi, sin_delta) * cos(degrees_to_radians(delta_alpha_prime_degrees)),
        mul_add(x * sin_xi, -cos_h, cos_delta),
    ));

    // 3.12. Calculate the topocentric local hour angle, H' (in degrees)
    let h_prime_degrees = h_degrees - delta_alpha_prime_degrees;
    let delta_prime = degrees_to_radians(delta_prime_degrees);
    let h_prime = degrees_to_radians(h_prime_degrees);
    let (sin_delta_prime, cos_delta_prime) = sin_cos(delta_prime);
    let (sin_h_prime, cos_h_prime) = sin_cos(h_prime);

    // 3.13. Calculate the topocentric zenith and azimuth angles
    let cos_zenith = mul_add(
        sin_phi,
        sin_delta_prime,
        cos_phi * cos_delta_prime * cos_h_prime,
    );
    // Roundoff can put the cosine just outside [-1, 1] at zenith or nadir.
    let zenith_angle = radians_to_degrees(acos(cos_zenith.clamp(-1.0, 1.0)));

    // 3.14. Calculate the topocentric azimuth angle
    let azimuth = normalize_degrees_0_to_360(
        180.0
            + radians_to_degrees(atan2(
                sin_h_prime,
                mul_add(
                    sin_delta_prime / cos_delta_prime,
                    -cos_phi,
                    cos_h_prime * sin_phi,
                ),
            )),
    );

    // Apply atmospheric refraction if requested
    let elevation_angle = 90.0 - zenith_angle;
    let final_zenith = refraction.map_or(zenith_angle, |correction| {
        if elevation_angle > -0.83337 {
            let pressure = correction.pressure();
            let temperature = correction.temperature();
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

#[allow(clippy::suboptimal_flops)] // Keep the vector geometry directly readable.
#[allow(clippy::unnecessary_wraps)] // Matches the fallible position-provider contract.
pub fn event_position(time: JulianDate, location: Location) -> Result<EventPosition> {
    let Location {
        latitude,
        longitude,
    } = location;
    let phi = latitude.to_radians();
    let u = atan(EARTH_FLATTENING_FACTOR * tan(phi));
    let (sin_phi, cos_phi) = sin_cos(phi);
    let x = cos(u);
    let z = EARTH_FLATTENING_FACTOR * sin(u);
    let parts = time_dependent(time);
    let delta = parts.delta_degrees.to_radians();
    let hour_angle = parts.nu_degrees + longitude - parts.alpha_degrees;
    let parallax = sin((8.794 / (3600.0 * parts.r)).to_radians());
    // SPA's topocentric parallax correction in vector form, stable at the zenith.
    let vx = cos(delta) * cos(hour_angle.to_radians()) - x * parallax;
    let vy = cos(delta) * sin(hour_angle.to_radians());
    let vz = sin(delta) - z * parallax;
    let projection = (cos_phi * vx + sin_phi * vz) / sqrt(vx * vx + vy * vy + vz * vz);
    Ok(EventPosition {
        elevation: asin(projection.clamp(-1.0, 1.0)).to_degrees(),
        hour_angle,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn angular_distance(a: f64, b: f64) -> f64 {
        let diff = (a - b).abs();
        diff.min(360.0 - diff)
    }

    #[test]
    fn test_time_dependent_tracks_seasonal_geometry() {
        let june_solstice =
            time_dependent(JulianDate::from_utc(2023, 6, 21, 12, 0, 0.0, 69.0).unwrap());
        let december_solstice =
            time_dependent(JulianDate::from_utc(2023, 12, 22, 12, 0, 0.0, 69.0).unwrap());
        let march_equinox =
            time_dependent(JulianDate::from_utc(2023, 3, 20, 12, 0, 0.0, 69.0).unwrap());

        assert!(june_solstice.delta_degrees > 23.0);
        assert!(june_solstice.delta_degrees < 24.0);
        assert!(angular_distance(june_solstice.alpha_degrees, 90.0) < 2.0);

        assert!(december_solstice.delta_degrees < -23.0);
        assert!(december_solstice.delta_degrees > -24.0);
        assert!(angular_distance(december_solstice.alpha_degrees, 270.0) < 2.0);

        assert!(march_equinox.delta_degrees.abs() < 1.0);
    }

    #[test]
    fn test_time_dependent_earth_radius_vector_changes_over_year() {
        let near_perihelion =
            time_dependent(JulianDate::from_utc(2023, 1, 4, 12, 0, 0.0, 69.0).unwrap());
        let near_aphelion =
            time_dependent(JulianDate::from_utc(2023, 7, 4, 12, 0, 0.0, 69.0).unwrap());

        assert!(near_perihelion.r > 0.98);
        assert!(near_perihelion.r < 0.99);
        assert!(near_aphelion.r > 1.01);
        assert!(near_aphelion.r < 1.02);
        assert!(near_aphelion.r > near_perihelion.r);
    }
}
