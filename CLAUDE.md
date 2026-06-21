# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

QuantumControlTestUtils.jl collects methods used for testing and benchmarking within the JuliaQuantumControl organization. It provides a test runner with coverage reporting, generators for random quantum objects, a dummy optimization problem/method for exercising the optimization machinery, and a logger for capturing log messages in tests.

It is part of the JuliaQuantumControl ecosystem and depends on QuantumControl.jl.

## Development Commands

### Development Environment
- `make devrepl`: Start an interactive REPL for testing and building documentation (recommended)
- `julia -i --project=test devrepl.jl`: Alternative way to start the development REPL

### Testing
- `make test`: Run the complete test suite
- `julia --project=test --banner=no --startup-file=yes -e 'include("devrepl.jl"); test()'`: Alternative test command

### Documentation
- `make docs`: Build the documentation

### Code Quality
- `make codestyle`: Apply JuliaFormatter to the entire project
- Requires `../.JuliaFormatter.toml` (org-level configuration)
- The JuliaFormatter version is pinned (`=2.3.0` in `test/Project.toml`); CI rejects formatting from other versions

### Cleanup
- `make clean`: Clean up build/doc/testing artifacts
- `make distclean`: Restore to a clean checkout state

## Architecture

### Module Structure
- `src/QuantumControlTestUtils.jl`: Main module; includes the components below
- `src/runner.jl`: Test runner (`test`) and coverage helpers (`show_coverage`, `collect_coverage`, `generate_coverage_html`)
- `src/random.jl`: `RandomObjects` submodule with generators for random states, matrices, and dynamic generators
- `src/dummy_problem.jl`: `DummyOptimization` submodule providing `dummy_control_problem` and `optimize_with_dummy_method`
- `src/logging.jl`: `QuantumTestLogger` for capturing and asserting on log messages in tests

### Key Components

#### Test Runner (`test`)
Runs a package test suite in a subprocess with coverage collection. Used by the `make test` target across the organization's packages.

#### Random Objects (`RandomObjects`)
Generators such as `random_state_vector`, `random_matrix`, and `random_dynamic_generator`, plus dense/sparse and Hermitian/real variants. All generators accept an `rng` keyword for reproducibility; prefer passing an explicit `rng` (e.g. a `StableRNG`) in tests rather than relying on the global RNG.

#### Dummy Optimization (`DummyOptimization`)
`dummy_control_problem` builds a `ControlProblem` with random components, and `optimize_with_dummy_method` runs a trivial optimization. These exercise the `QuantumControl.optimize` machinery without a real optimization method.

## Changelog

`CHANGELOG.md` follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) / [SemVer](https://semver.org/). Non-obvious conventions:

* Record user-facing changes under `## [Unreleased]` as bullets with an inline category prefix (`Added:`/`Changed:`/`Deprecated:`/`Removed:`/`Fixed:`/`Security:`), not `###` subsections; link issues/PRs as `[[#123]]`, issue before its resolving PR (`[[#91], [#93]]`). Exclude CI, dependency bumps, formatting, and internal-only changes — a leading underscore (e.g. `_helper`) marks a name as internal.
* Pre-1.0, Julia treats every `v0.x.0` as breaking, so non-breaking changes go into a `v0.x.y` bugfix release.
* Version links point to the release page (`[vX.Y.Z]: …/releases/tag/vX.Y.Z`); only `[Unreleased]` uses a compare link (`…/compare/v<latest>..HEAD`).
* `pull/` vs `issues/` can't be verified by loading the URL (GitHub redirects between them); confirm the category with `gh api repos/JuliaQuantumControl/QuantumControlTestUtils.jl/issues/<N> --jq 'if has("pull_request") then "pull" else "issue" end'`.
* Releasing on a `release-*` branch: rename `## [Unreleased]` to `## [vX.Y.Z] — YYYY-MM-DD` and point `[Unreleased]` at `…/compare/vX.Y.Z..HEAD`, but do **not** add a fresh `## [Unreleased]` heading — re-add it when merging back to `master`.
* `make check-changelog` validates links (textual, no network; also run in CI via `make codestyle`); `make changelog` additionally fills in missing `[#N]` targets, so you can just write `[[#123]]`. Neither verifies that links resolve — check that, and the issue/PR category, manually.

## Development Notes

- Part of the JuliaQuantumControl ecosystem; uses shared development scripts in `../scripts/` (notably `installorg.jl`)
- Designed to work within the JuliaQuantumControl development environment
- Tests automatically use the current dev versions of sibling packages
- Code formatting follows JuliaQuantumControl organization standards
- Uses `devrepl.jl` for unified development environment setup
