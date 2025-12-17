# Changelog

## Unreleased

### Fixed

- Fix missing sunrise/twilight events near midnight at timezone boundaries in `chrono` sunrise/sunset APIs.
- Keep sunrise/sunset bracketing transit on the requested local date even when SPA wraps near UTC midnight.

## [0.4.4] - 2025-11-30

### Changed

- Trim `SpaTimeDependent` to only keep values used by `spa_with_time_dependent_parts()`.
- Align single and multi sunrise behaviour in polar conditions. This may result in output changes.

## [0.4.3] - 2025-11-14

### Changed

- fix: make `HoursUtc::day_and_hours` normalise multi-day offsets so the hour component stays within `[0, 24)`.

## [0.4.2] - 2025-10-26

### Fixed

- Prevent local timezone inputs from jumping to the previous civil day when truncating to midnight in SPA calculations.

## [0.4.1] - 2025-10-25

### Added

- More input validation.

### Changed

- Adapt delta T behaviour in sunrise/sunset calculation to match SPA paper more closely.

## [0.4.0] - 2025-10-07

### Added

- Basic sunrise/sunset/twilight calculations without chrono

### Changed

- **BREAKING** Make API more consistent (re-exports changed)

## [0.3.6] - 2025-10-04

### Changed

- reduce dependency footprint and made `chrono` optional
- drop down to Rust 2021 edition (MSRV 1.70)

## [0.3.5] - 2025-10-03

### Changed

- build, documentation, and test fixes

## [0.3.4] - 2025-09-24

### Added

- feat: basic no_std support

## [0.3.3] - 2025-09-23

### Changed

- **BREAKING** fix: remove useless parameters from `spa_with_time_dependent_parts()` and make `SpaTimeDependent` opaque

## [0.3.1] - 2025-09-22

### Changed

- **BREAKING** fix: `DeltaT::estimate_from_date_like()` now takes parameter by value instead of reference

## [0.3.0] - 2025-09-22

### Added

- `RefractionCorrection` struct with validation for atmospheric correction parameters
- Split SPA functions `spa_time_dependent_parts()` and `spa_with_time_dependent_parts()` for performance optimization of coordinate sweeps (same datetime, multiple coordinates)

### Changed

- **BREAKING**: `solar_position()` refraction parameter changed from separate pressure/temperature to `Option<RefractionCorrection>`
- **BREAKING**: Removed `solar_position_no_refraction()` function
- **BREAKING**: `grena3::solar_position()` refraction parameter unified with SPA API

## [0.2.2] - 2025-09-18

### Added

- `DeltaT::estimate_from_date_like()` - convenience method accepting any chrono `Datelike` type

### Fixed

- Polar latitude transit accuracy - sunrise/sunset calculations now apply full iterative corrections for polar conditions, fixing 14-second error at extreme latitudes like Svalbard

## [0.2.1] - 2025-09-16

### Added

- Initial port from Java version
- NREL SPA algorithm
- Grena3 algorithm
- Sunrise/sunset/twilight calculations
- Support for all twilight types
- Delta T estimation
- Atmospheric refraction corrections
- Comprehensive test coverage
- Thread-safe API
