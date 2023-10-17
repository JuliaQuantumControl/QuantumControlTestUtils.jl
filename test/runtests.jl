using Test
using SafeTestsets

@time @testset verbose = true "QuantumControlTestUtils" begin

    println("\n* Random objects (test_random.jl):")
    @time @safetestset "Random objects" begin
        include("test_random.jl")
    end

    println("\n* Doctests (test_doctest.jl):")
    @time @safetestset "Doctests" begin
        include("test_doctest.jl")
    end

    print("\n")

end
