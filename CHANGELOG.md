# Release Notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [v0.3.2] — 2026-06-21

* Fixed: The random-object generators now correctly use the explicitly provided `rng` when generating random pulses
* Fixed: `random_hermitian_real_matrix` now correctly uses the explicitly provided `rng`

## [v0.3.1] — 2024-09-04

* Changed: Depend directly on `QuantumControl` instead of `QuantumControlBase` [[#23]]

## [v0.3.0] — 2024-01-23

* Changed: Adapted to the renaming of "objective" to "trajectory" (`Objective` → `Trajectory`) in `QuantumControl`
* Changed: Depend on `QuantumControlBase` [[#21]]

## [v0.2.1] — 2023-12-02

* Changed: `random_dynamic_generator` now allows custom amplitudes
* Changed: `dummy_control_problem` now explicitly specifies a `prop_method`
* Changed: The `test` runner now also collects coverage data from extension modules

## [v0.2.0] — 2023-10-18

* Changed: The minimum supported Julia version is now 1.9
* Changed: The `test` runner now collects coverage data in the `.coverage` directory

## [v0.1.5] — 2023-05-16

* Changed: Documentation and test improvements

## [v0.1.4] — 2023-04-04

* Changed: Compatibility with `QuantumControl` 0.7

## [v0.1.3] — 2023-03-19

* Added: `dummy_control_problem` now accepts a `pulses_as_controls` keyword argument
* Changed: `dummy_control_problem` is now reproducible

## [v0.1.2] — 2023-02-16

* Changed: Compatibility with `QuantumControl` 0.6

## [v0.1.1] — 2023-02-15

* Added: `DummyOptimization` submodule, providing `dummy_control_problem` and `optimize_with_dummy_method`
* Added: `random_dynamic_generator`
* Changed: Random matrices now use an exact spectral radius

## [v0.1.0] — 2022-12-21

Initial public release, providing the `test` runner, the `RandomObjects` submodule (`random_state_vector`, `random_matrix`, and variants), and the `QuantumTestLogger`.

[Unreleased]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/compare/v0.3.2..HEAD
[v0.3.2]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.3.2
[v0.3.1]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.3.1
[v0.3.0]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.3.0
[v0.2.1]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.2.1
[v0.2.0]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.2.0
[v0.1.5]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.1.5
[v0.1.4]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.1.4
[v0.1.3]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.1.3
[v0.1.2]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.1.2
[v0.1.1]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.1.1
[v0.1.0]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/releases/tag/v0.1.0
[#21]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/pull/21
[#23]: https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/pull/23
