using QuantumControlTestUtils
using Documenter

DocMeta.setdocmeta!(
    QuantumControlTestUtils,
    :DocTestSetup,
    :(using QuantumControlTestUtils);
    recursive=true
)

PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl"

println("Starting makedocs")

makedocs(;
    authors=AUTHORS,
    sitename="QuantumControlTestUtils.jl",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://JuliaQuantumControl.github.io/QuantumControlTestUtils.jl",
        assets=String[],
        footer="[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)."
    ),
    pages=["Home" => "index.md",]
)

println("Finished makedocs")

deploydocs(; repo="github.com/JuliaQuantumControl/QuantumControlTestUtils.jl")
