using QuantumControlTestUtils
using Documenter

DocMeta.setdocmeta!(
    QuantumControlTestUtils,
    :DocTestSetup,
    :(using QuantumControlTestUtils);
    recursive=true
)

makedocs(;
    modules=[QuantumControlTestUtils],
    authors="Michael Goerz <mail@michaelgoerz.net",
    repo="https://github.com/JuliaQuantumControl/QuantumControlTestUtils.jl/blob/{commit}{path}#{line}",
    sitename="QuantumControlTestUtils.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaQuantumControl.github.io/QuantumControlTestUtils.jl",
        edit_link="master",
        assets=String[]
    ),
    pages=["Home" => "index.md",]
)

deploydocs(;
    repo="github.com/JuliaQuantumControl/QuantumControlTestUtils.jl",
    devbranch="master"
)
