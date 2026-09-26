//! Solar positions and reusable calculations for a fixed time.

use crate::{
    Error, Location, RefractionCorrection, Result, SolarPosition, grena3, spa, time::JulianDate,
};
#[cfg(feature = "chrono")]
use chrono::{DateTime, TimeZone};

/// A reusable solar-position calculator, using SPA by default.
///
/// Height is in metres above sea level; use zero for sea level. Grena3 requires
/// zero height. Pass `None` for an unrefracted position, or `Some` atmospheric
/// conditions to apply the selected model's refraction correction.
///
/// SPA is designed for years -2000 through 6000, and Grena3 for 2010 through 2110.
/// Instances are immutable, `Send` and `Sync`, and require no allocation.
#[derive(Debug, Clone, Copy, Default)]
pub struct SolarPositions {
    model: Model,
}

#[derive(Debug, Clone, Copy, Default)]
enum Model {
    #[default]
    Spa,
    Grena3,
}

impl SolarPositions {
    /// Creates a calculator using SPA.
    #[must_use]
    pub const fn new() -> Self {
        Self { model: Model::Spa }
    }

    /// Creates a calculator using Grena3.
    #[must_use]
    pub const fn grena3() -> Self {
        Self {
            model: Model::Grena3,
        }
    }

    /// Calculates a position at a timezone-aware timestamp.
    ///
    /// `height` is in metres and `delta_t` is TT minus UT1 in seconds.
    /// UTC approximates UT1. The timestamp is borrowed, so any chrono time zone works
    /// without requiring the caller to copy or clone it.
    ///
    /// # Errors
    /// Returns an error for an invalid timestamp, non-finite delta T or height,
    /// invalid coordinates, nonzero height with Grena3, or a non-finite result.
    #[cfg(feature = "chrono")]
    #[cfg_attr(docsrs, doc(cfg(feature = "chrono")))]
    #[inline]
    pub fn at<Tz: TimeZone>(
        &self,
        time: &DateTime<Tz>,
        location: Location,
        height: f64,
        delta_t: f64,
        refraction: Option<RefractionCorrection>,
    ) -> Result<SolarPosition> {
        self.for_time(time, delta_t)?
            .at(location, height, refraction)
    }

    /// Calculates a position from continuous UT Julian time, which includes delta T.
    ///
    /// `height` is in metres. This method needs neither chrono nor allocation.
    ///
    /// # Errors
    /// Returns an error for invalid coordinates, non-finite height, nonzero height
    /// with Grena3, or a non-finite result.
    #[inline]
    pub fn at_from_julian(
        &self,
        time: JulianDate,
        location: Location,
        height: f64,
        refraction: Option<RefractionCorrection>,
    ) -> Result<SolarPosition> {
        self.for_time_from_julian(time)
            .at(location, height, refraction)
    }

    /// Prepares time-dependent calculations for reuse at multiple locations.
    ///
    /// `delta_t` is TT minus UT1 in seconds. The result owns its cached values and
    /// does not borrow the timestamp or calculator.
    ///
    /// # Errors
    /// Returns an error for an invalid timestamp or non-finite delta T.
    #[cfg(feature = "chrono")]
    #[cfg_attr(docsrs, doc(cfg(feature = "chrono")))]
    #[inline]
    pub fn for_time<Tz: TimeZone>(
        &self,
        time: &DateTime<Tz>,
        delta_t: f64,
    ) -> Result<PreparedPositions> {
        Ok(self.for_time_from_julian(JulianDate::from_datetime(time, delta_t)?))
    }

    /// Prepares calculations from continuous UT Julian time, which includes delta T.
    ///
    /// This method needs neither chrono nor allocation.
    #[must_use]
    #[inline]
    pub fn for_time_from_julian(&self, time: JulianDate) -> PreparedPositions {
        let model = match self.model {
            Model::Spa => PreparedModel::Spa(spa::time_dependent(time)),
            Model::Grena3 => PreparedModel::Grena3(grena3::time_dependent(time)),
        };
        PreparedPositions { model }
    }
}

/// Solar positions at one time, prepared by [`SolarPositions::for_time_from_julian`].
///
/// Both models reuse their time-dependent calculations. The cached values are owned,
/// so this value can outlive its calculator and input timestamp. It is immutable,
/// `Send` and `Sync`, and requires no allocation.
#[derive(Debug, Clone, Copy)]
pub struct PreparedPositions {
    model: PreparedModel,
}

#[derive(Debug, Clone, Copy)]
enum PreparedModel {
    Spa(spa::TimeDependent),
    Grena3(grena3::TimeDependent),
}

impl PreparedPositions {
    /// Calculates a position at this prepared time.
    ///
    /// `height` is in metres above sea level; use zero for sea level. Grena3 requires
    /// zero height. Pass `None` to omit refraction correction. Height and atmospheric
    /// conditions can vary between locations.
    ///
    /// # Errors
    /// Returns an error for invalid coordinates, non-finite height, nonzero height
    /// with Grena3, or a non-finite result.
    #[inline]
    pub fn at(
        &self,
        location: Location,
        height: f64,
        refraction: Option<RefractionCorrection>,
    ) -> Result<SolarPosition> {
        if !height.is_finite() {
            return Err(Error::InvalidHeight { value: height });
        }
        match &self.model {
            PreparedModel::Spa(parts) => spa::solar_position(
                location.latitude,
                location.longitude,
                height,
                refraction,
                parts,
            ),
            PreparedModel::Grena3(parts) => {
                if height != 0.0 {
                    return Err(Error::InvalidHeight { value: height });
                }
                grena3::solar_position(location.latitude, location.longitude, refraction, parts)
            }
        }
    }
}
