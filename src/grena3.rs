//! Grena3 algorithm implementation.
//!
//! Implements algorithm #3 from Grena (2012).
//!
//! Designed for 2010-2110 with 0.01° accuracy. ~10x faster than SPA.
//!
//! Reference: Grena, R. (2012). Five new algorithms for the computation of sun position from 2010 to 2110.
//! Solar Energy, 86(5), 1323-1337. DOI: <http://dx.doi.org/10.1016/j.solener.2012.01.024>

#![allow(clippy::unreadable_literal)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::suboptimal_flops)]

use crate::error::check_coordinates;
use crate::events::EventPosition;
use crate::math::{
    PI, asin, atan2, cos, degrees_to_radians, normalize_degrees_0_to_360, radians_to_degrees,
    rem_euclid, sin, sqrt, tan,
};
use crate::time::JulianDate;
use crate::{Location, RefractionCorrection, Result, SolarPosition};

/// Angles in radians, independent of observer location.
#[derive(Debug, Clone, Copy)]
pub struct TimeDependent {
    sidereal_time: f64,
    alpha: f64,
    sin_delta: f64,
    cos_delta: f64,
}

#[inline]
pub fn solar_position(
    latitude: f64,
    longitude: f64,
    refraction: Option<RefractionCorrection>,
    parts: &TimeDependent,
) -> Result<SolarPosition> {
    check_coordinates(latitude, longitude)?;

    let Position {
        elevation: e_p,
        azimuth: gamma,
        ..
    } = position(parts, latitude, longitude);

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

// Unrefracted angles in radians; azimuth is measured westwards from south.
struct Position {
    elevation: f64,
    azimuth: f64,
    hour_angle: f64,
}

#[inline]
pub fn time_dependent(time: JulianDate) -> TimeDependent {
    // Continuous UT days from Grena's epoch (2060-01-01), avoiding the calendar
    // formula's rounded hours-to-days constant and its tiny midnight jump.
    let t = time.julian_date() - 2_473_459.5;
    let t_e = t + 1.1574e-5 * time.delta_t();
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
    let alpha = rem_euclid(atan2(s_lambda * c_epsilon, c_lambda), 2.0 * PI);

    // Calculate declination (delta)
    let delta = asin(s_lambda * s_epsilon);

    let sin_delta = sin(delta);
    let cos_delta = sqrt(1.0 - sin_delta * sin_delta);
    TimeDependent {
        sidereal_time: 1.7528311 + 6.300388099 * t,
        alpha,
        sin_delta,
        cos_delta,
    }
}

fn position(parts: &TimeDependent, latitude: f64, longitude: f64) -> Position {
    // Calculate hour angle (H)
    let h = rem_euclid(
        parts.sidereal_time + degrees_to_radians(longitude) - parts.alpha + PI,
        2.0 * PI,
    ) - PI;

    // Calculate topocentric coordinates
    let s_phi = sin(degrees_to_radians(latitude));
    let c_phi = sqrt(1.0 - s_phi * s_phi);
    let s_delta = parts.sin_delta;
    let c_delta = parts.cos_delta;
    let s_h = sin(h);
    let c_h = cos(h);

    // Roundoff can put the sine just outside [-1, 1] at the zenith or nadir.
    let s_epsilon0 = (s_phi * s_delta + c_phi * c_delta * c_h).clamp(-1.0, 1.0);
    let e_p = asin(s_epsilon0) - 4.26e-5 * sqrt(1.0 - s_epsilon0 * s_epsilon0);
    let gamma = atan2(s_h, c_h * s_phi - s_delta * c_phi / c_delta);

    Position {
        elevation: e_p,
        azimuth: gamma,
        hour_angle: h,
    }
}

#[allow(clippy::unnecessary_wraps)] // Matches the fallible position-provider contract.
pub fn event_position(time: JulianDate, location: Location) -> Result<EventPosition> {
    let Location {
        latitude,
        longitude,
    } = location;
    let position = position(&time_dependent(time), latitude, longitude);
    Ok(EventPosition {
        elevation: position.elevation.to_degrees(),
        hour_angle: position.hour_angle.to_degrees(),
    })
}
