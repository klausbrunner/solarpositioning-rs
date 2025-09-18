# Changelog

## [Unreleased]

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
