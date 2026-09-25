use chrono::{DateTime, Utc};
use solar_positioning::{spa, SunriseResult};

// The reference tables list events within a UTC calendar day. The library pairs
// events around transit, so a table's rise or set may belong to an adjacent transit.
pub fn events_on_utc_date(
    date: DateTime<Utc>,
    latitude: f64,
    longitude: f64,
    elevation: f64,
) -> solar_positioning::Result<SunriseResult<DateTime<Utc>>> {
    let calculate = |date| spa::sunrise_sunset(date, latitude, longitude, 0.0, elevation);
    let mut result = calculate(date)?;
    if let SunriseResult::RegularDay {
        sunrise, sunset, ..
    } = &mut result
    {
        for (event, rising) in [(sunrise, true), (sunset, false)] {
            let offset = date.date_naive() - event.date_naive();
            if !offset.is_zero() {
                let adjacent = calculate(date + offset)?;
                *event = *if rising {
                    adjacent.sunrise()
                } else {
                    adjacent.sunset()
                }
                .ok_or(solar_positioning::Error::ComputationError {
                    message:
                        "SPA classifies the adjacent transit as polar; reference event unavailable",
                })?;
            }
            assert_eq!(
                event.date_naive(),
                date.date_naive(),
                "wrong reference event date"
            );
        }
    }
    Ok(result)
}
