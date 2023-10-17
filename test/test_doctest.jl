using Test
using Documenter: DocMeta, doctest
using QuantumControlTestUtils

@testset "run doctest" begin
    DocMeta.setdocmeta!(
        QuantumControlTestUtils,
        :DocTestSetup,
        :(using QuantumControlTestUtils);
        recursive=true,
        warn=false
    )
    doctest(QuantumControlTestUtils)
end
