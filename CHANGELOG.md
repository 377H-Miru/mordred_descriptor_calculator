# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-06-29

### Added
- Professional CLI entry point `mordred-desc`.
- Package structure using standard `src` layout (`mordred_descriptor_calculator`).
- Support for CLI arguments replacing mandatory config file requirement.
- Multi-stage failure tracking and structured error logging (`output.csv.errors.csv`).
- Configurable parallel worker process management and safe sequential fallback.
- Support for custom conjugation descriptors and 2D/3D computation modes (`--only-2d`, `--include-3d`, `--include-conjugation`).
- Full test suite with unit tests and CLI integration tests.
- Added GitHub Actions CI workflow for Python 3.10.
