//! Appendix A.2 conformance, not a comparison against observed sunrise.
//!
//! Values generated with pvlib 0.15.2, transit_sunrise_sunset, delta_ut1 = 0,
//! horizon = -0.8333 degrees. The large delta_t case exposes a duplicated TT offset.
//! https://github.com/pvlib/pvlib-python/blob/v0.15.2/pvlib/spa.py
//! Cases avoid day wrapping, which is tested separately.
//! The polar case deliberately retains the paper's single-correction approximation.

use solar_positioning::{spa, SunriseResult};

#[test]
fn sunrise_interpolation_uses_tt_midnight() {
    // date, latitude, longitude, delta_t, [transit, sunrise, sunset] in UTC hours
    let cases = [
        (
            (2003, 10, 17),
            39.742476,
            0.0,
            67.0,
            [11.757134843700461, 6.195069657795959, 17.309518482850656],
        ),
        (
            (2024, 3, 20),
            48.21,
            16.37,
            69.184,
            [11.030678879088825, 4.94628897620572, 17.13000080221229],
        ),
        (
            (2024, 8, 30),
            80.0,
            0.0,
            69.184,
            [12.007834764719009, 0.9303217564026515, 22.527126493785115],
        ),
        (
            (2024, 3, 20),
            48.21,
            16.37,
            3600.0,
            [11.033161676393615, 4.947567514975866, 17.13368758128749],
        ),
    ];
    for ((year, month, day), latitude, longitude, delta_t, expected) in cases {
        let SunriseResult::RegularDay {
            sunrise,
            transit,
            sunset,
        } = spa::sunrise_sunset_utc(year, month, day, latitude, longitude, delta_t, -0.8333)
            .unwrap()
        else {
            panic!("expected regular day");
        };
        for (actual, expected) in [transit, sunrise, sunset].into_iter().zip(expected) {
            let error_seconds = (actual.hours() - expected).abs() * 3600.0;
            assert!(
                error_seconds < 0.001,
                "{year}-{month}-{day}, delta_t={delta_t}: error {error_seconds}s"
            );
        }
    }
}
