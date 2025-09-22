# Changelog

## [Unreleased]

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
