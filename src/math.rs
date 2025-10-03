//! Mathematical utilities for solar position calculations.

#![allow(clippy::many_single_char_names)]

#[cfg(not(feature = "std"))]
use libm;

/// Mathematical constants
pub const PI: f64 = core::f64::consts::PI;

/// Converts degrees to radians.
#[inline]
pub const fn degrees_to_radians(degrees: f64) -> f64 {
    degrees.to_radians()
}

/// Converts radians to degrees.
#[inline]
pub const fn radians_to_degrees(radians: f64) -> f64 {
    radians.to_degrees()
}

/// Normalizes an angle in degrees to the range [0, 360).
pub fn normalize_degrees_0_to_360(degrees: f64) -> f64 {
    let normalized = degrees % 360.0;
    if normalized < 0.0 {
        normalized + 360.0
    } else {
        normalized
    }
}

/// Computes a polynomial using Horner's method for numerical stability.
///
/// Coefficients are ordered [a₀, a₁, a₂, ...] for a₀ + a₁x + a₂x² + ...
pub fn polynomial(coeffs: &[f64], x: f64) -> f64 {
    let Some(&last) = coeffs.last() else {
        return 0.0;
    };

    // Horner's method: reverse iteration for numerical stability
    let mut result = last;
    for &coeff in coeffs.iter().rev().skip(1) {
        result = mul_add(result, x, coeff);
    }
    result
}

/// Computes sin(x) using the appropriate function for the compilation target.
#[inline]
pub fn sin(x: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.sin();

    #[cfg(not(feature = "std"))]
    return libm::sin(x);
}

/// Computes cos(x) using the appropriate function for the compilation target.
#[inline]
pub fn cos(x: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.cos();

    #[cfg(not(feature = "std"))]
    return libm::cos(x);
}

/// Computes tan(x) using the appropriate function for the compilation target.
#[inline]
pub fn tan(x: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.tan();

    #[cfg(not(feature = "std"))]
    return libm::tan(x);
}

/// Computes asin(x) using the appropriate function for the compilation target.
#[inline]
pub fn asin(x: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.asin();

    #[cfg(not(feature = "std"))]
    return libm::asin(x);
}

/// Computes acos(x) using the appropriate function for the compilation target.
#[inline]
pub fn acos(x: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.acos();

    #[cfg(not(feature = "std"))]
    return libm::acos(x);
}

/// Computes atan(x) using the appropriate function for the compilation target.
#[inline]
pub fn atan(x: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.atan();

    #[cfg(not(feature = "std"))]
    return libm::atan(x);
}

/// Computes atan2(y, x) using the appropriate function for the compilation target.
#[inline]
pub fn atan2(y: f64, x: f64) -> f64 {
    #[cfg(feature = "std")]
    return y.atan2(x);

    #[cfg(not(feature = "std"))]
    return libm::atan2(y, x);
}

/// Computes sqrt(x) using the appropriate function for the compilation target.
#[inline]
pub fn sqrt(x: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.sqrt();

    #[cfg(not(feature = "std"))]
    return libm::sqrt(x);
}

/// Computes floor(x) using the appropriate function for the compilation target.
#[inline]
pub fn floor(x: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.floor();

    #[cfg(not(feature = "std"))]
    return libm::floor(x);
}

/// Computes (x * a) + b with only one rounding error (fused multiply-add).
#[inline]
pub fn mul_add(x: f64, a: f64, b: f64) -> f64 {
    #[cfg(feature = "std")]
    return x.mul_add(a, b);

    #[cfg(not(feature = "std"))]
    return libm::fma(x, a, b);
}

/// Computes x^n for integer n.
#[inline]
pub fn powi(x: f64, n: i32) -> f64 {
    #[cfg(feature = "std")]
    return x.powi(n);

    #[cfg(not(feature = "std"))]
    return libm::pow(x, f64::from(n));
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_degree_radian_conversion() {
        assert!((degrees_to_radians(180.0) - PI).abs() < EPSILON);
        assert!((degrees_to_radians(90.0) - PI / 2.0).abs() < EPSILON);
        assert!((degrees_to_radians(0.0)).abs() < EPSILON);

        assert!((radians_to_degrees(PI) - 180.0).abs() < EPSILON);
        assert!((radians_to_degrees(PI / 2.0) - 90.0).abs() < EPSILON);
        assert!((radians_to_degrees(0.0)).abs() < EPSILON);
    }

    #[test]
    fn test_normalize_degrees_0_to_360() {
        assert_eq!(normalize_degrees_0_to_360(0.0), 0.0);
        assert_eq!(normalize_degrees_0_to_360(90.0), 90.0);
        assert_eq!(normalize_degrees_0_to_360(360.0), 0.0);
        assert_eq!(normalize_degrees_0_to_360(450.0), 90.0);
        assert_eq!(normalize_degrees_0_to_360(-90.0), 270.0);
        assert_eq!(normalize_degrees_0_to_360(-360.0), 0.0);
    }

    #[test]
    fn test_polynomial() {
        // Test empty coefficients
        assert_eq!(polynomial(&[], 5.0), 0.0);

        // Test constant polynomial
        assert_eq!(polynomial(&[3.0], 5.0), 3.0);

        // Test linear polynomial: 2 + 3x
        assert_eq!(polynomial(&[2.0, 3.0], 4.0), 14.0);

        // Test quadratic polynomial: 1 + 2x + 3x²
        assert!((polynomial(&[1.0, 2.0, 3.0], 2.0) - 17.0).abs() < EPSILON);
    }

    #[test]
    fn test_trigonometric_functions() {
        // Basic smoke tests - the actual implementation will depend on features
        assert!((sin(0.0)).abs() < EPSILON);
        assert!((cos(0.0) - 1.0).abs() < EPSILON);
        assert!((tan(0.0)).abs() < EPSILON);
    }
}
