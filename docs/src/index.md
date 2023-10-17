# QuantumControlTestUtils

```@eval
using Markdown
using Pkg

VERSION = Pkg.dependencies()[Base.UUID("d3fd27c9-1dfb-4e67-b0c0-90d0d87a1e48")].version

github_badge = "[![Github](https://img.shields.io/badge/JuliaQuantumControl-QuantumControlTestUtils.jl-blue.svg?logo=github)](https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl)"

version_badge = "![v$VERSION](https://img.shields.io/badge/version-v$VERSION-green.svg)"

Markdown.parse("$github_badge $version_badge")
```

The [QuantumControlTestUtils](https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl) package collects methods that are used for testing and benchmarking within the [JuliaQuantumControl](https://github.com/JuliaQuantumControl) organization

```@autodocs
Modules = [QuantumControlTestUtils, QuantumControlTestUtils.RandomObjects, QuantumControlTestUtils.DummyOptimization]
Private = false
```
