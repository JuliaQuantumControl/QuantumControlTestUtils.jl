# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

QuantumControlTestUtils.jl provides generators for random quantum objects (states, matrices, dynamic generators) used for testing and benchmarking within the JuliaQuantumControl organization.

The package deliberately depends only on the Julia standard library (LinearAlgebra, Random, SparseArrays). It appears in the test environments of many org packages, so any other dependency could restrict package versions in those environments. Do not add dependencies. In particular, it must not depend on `QuantumPropagators` or `QuantumControl`: `random_dynamic_generator` returns the terms of a generator, which callers pass to `hamiltonian`.

## Development Commands

Run `make help` for all targets. The development workflow is documented in the org-wide [CONTRIBUTING.md](https://github.com/JuliaQuantumControl/.github/blob/master/CONTRIBUTING.md) (`../.github/CONTRIBUTING.md` in the development environment).

- `make test`: Run the test suite in the `test` environment
- `make devrepl`: REPL with the `test` environment active and the `docs` environment stacked; run `include("test/runtests.jl")` or `include("docs/make.jl")`
- `make docs`: Build the documentation in the `docs` environment
- `make coverage` / `make htmlcoverage`: Test coverage
- `make codestyle`: Apply JuliaFormatter (version pinned in the `Makefile`) and check `CHANGELOG.md` and `[sources]`
- `make clean` / `make distclean`

The `test` and `docs` environments reference the package itself via `[sources]` (`{path = ".."}`); this needs Julia ≥ 1.11.

## Architecture

- `src/QuantumControlTestUtils.jl`: Main module
- `src/random.jl`: `RandomObjects` submodule with `random_state_vector`, `random_matrix`, and `random_dynamic_generator`, plus dense/sparse and Hermitian/real variants. All generators accept an `rng` keyword for reproducibility; prefer passing an explicit `rng` (e.g. a `StableRNG`) in tests rather than relying on the global RNG.

## Changelog

`CHANGELOG.md` follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) / [SemVer](https://semver.org/). Non-obvious conventions:

* Record user-facing changes under `## [Unreleased]` as bullets with an inline category prefix (`Added:`/`Changed:`/`Deprecated:`/`Removed:`/`Fixed:`/`Security:`), not `###` subsections; link issues/PRs as `[[#123]]`, issue before its resolving PR (`[[#91], [#93]]`). Exclude CI, dependency bumps, formatting, and internal-only changes — a leading underscore (e.g. `_helper`) marks a name as internal.
* Pre-1.0, Julia treats every `v0.x.0` as breaking, so non-breaking changes go into a `v0.x.y` bugfix release.
* Version links point to the release page (`[vX.Y.Z]: …/releases/tag/vX.Y.Z`); only `[Unreleased]` uses a compare link (`…/compare/v<latest>..HEAD`).
* `pull/` vs `issues/` can't be verified by loading the URL (GitHub redirects between them); confirm the category with `gh api repos/JuliaQuantumControl/QuantumControlTestUtils.jl/issues/<N> --jq 'if has("pull_request") then "pull" else "issue" end'`.
* Releasing on a `release-*` branch: rename `## [Unreleased]` to `## [vX.Y.Z] — YYYY-MM-DD` and point `[Unreleased]` at `…/compare/vX.Y.Z..HEAD`, but do **not** add a fresh `## [Unreleased]` heading — re-add it when merging back to `master`.
* Separate sections (`## [...]` headings) with two blank lines. When merging a release back to `master`, the re-added `## [Unreleased]` section is empty, so its heading has two blank lines before and after it.
* `make check-changelog` validates links (textual, no network; also run in CI via `make codestyle`); `make changelog` additionally fills in missing `[#N]` targets, so you can just write `[[#123]]`. Neither verifies that links resolve — check that, and the issue/PR category, manually.
