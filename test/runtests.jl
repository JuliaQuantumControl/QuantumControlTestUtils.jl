using Test
using SafeTestsets

@time @testset verbose = true "QuantumControlTestUtils" begin

    print("\n* Random objects (test_random.jl):")
    @time @safetestset "Random objects" begin
        include("test_random.jl")
    end

    print("\n")

end
