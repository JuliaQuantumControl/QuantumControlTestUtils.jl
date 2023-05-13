using Test
using StableRNGs: StableRNG
using LinearAlgebra: norm, eigvals
import SparseArrays

using QuantumControl.Controls: evaluate
using QuantumControl.Generators: Generator
using QuantumControl.Interfaces: check_generator

using QuantumControlTestUtils.RandomObjects
# random_matrix, random_state_vector, random_dynamic_generator

N = 100

relerr(a, b) = abs(a - b) / max(abs(a), abs(b))

# How close do eigenvalues have to be to spectral radius? The larer N, the
# smaller rtol can be
≊(a, b) = isapprox(a, b; rtol=0.1)  # for larger N, rtol could be smaller

@testset "Random Hermitian Dense Real Matrix" begin

    H = random_matrix(N; hermitian=true, complex=false)
    @test H isa Matrix{Float64}
    @test size(H) == (N, N)

    rng = StableRNG(3405091510)

    ρ = 1.0
    H = random_matrix(N; hermitian=true, complex=false, rng)
    @test H isa Matrix{Float64}
    @test norm(H - H') ≈ 0.0
    @test norm(H - transpose(H)) ≈ 0.0
    λ = eigvals(H)
    @test λ isa Vector{Float64}
    sort!(λ)
    #@show relerr(λ[1], -ρ), relerr(λ[end], ρ)
    @test λ[1] ≊ -ρ
    @test λ[end] ≊ ρ

    ρ = 3.0
    H = random_matrix(N; hermitian=true, complex=false, spectral_radius=ρ, rng)
    @test H isa Matrix{Float64}
    λ = sort(eigvals(H))
    #@show relerr(λ[1], -ρ), relerr(λ[end], ρ)
    @test λ[1] ≊ -ρ
    @test λ[end] ≊ ρ

    ρ = 3.0
    H = random_matrix(
        N;
        hermitian=true,
        complex=false,
        spectral_radius=ρ,
        rng,
        exact_spectral_radius=true
    )
    @test H isa Matrix{Float64}
    λ = sort(eigvals(H))
    @test λ[1] ≈ -ρ
    @test λ[end] ≈ ρ

end


@testset "Random Hermitian Dense Complex Matrix" begin

    H = random_matrix(N; hermitian=true)
    @test H isa Matrix{ComplexF64}
    @test size(H) == (N, N)

    rng = StableRNG(3727996169)

    ρ = 1.0
    H = random_matrix(N; hermitian=true, rng)
    @test H isa Matrix{ComplexF64}
    @test norm(H - H') ≈ 0.0
    @test norm(H - transpose(H)) > 0.0
    λ = eigvals(H)
    @test λ isa Vector{Float64}
    sort!(λ)
    #@show relerr(λ[1], -ρ), relerr(λ[end], ρ)
    @test λ[1] ≊ -ρ
    @test λ[end] ≊ ρ

    ρ = 3.0
    H = random_matrix(N; hermitian=true, spectral_radius=ρ, rng)
    @test H isa Matrix{ComplexF64}
    λ = sort(eigvals(H))
    #@show relerr(λ[1], -ρ), relerr(λ[end], ρ)
    @test λ[1] ≊ -ρ
    @test λ[end] ≊ ρ

    ρ = 3.0
    H = random_matrix(N; hermitian=true, spectral_radius=ρ, rng, exact_spectral_radius=true)
    @test H isa Matrix{ComplexF64}
    λ = sort(eigvals(H))
    @test λ[1] ≈ -ρ
    @test λ[end] ≈ ρ

end


@testset "Random Hermitian Sparse Real Matrix" begin

    H = random_matrix(N; density=0.5, hermitian=true, complex=false)
    @test H isa SparseArrays.SparseMatrixCSC{Float64,Int64}
    @test size(H) == (N, N)

    rng = StableRNG(1572573603)

    ρ = 1.0
    d = 0.5
    H = random_matrix(N; density=d, hermitian=true, complex=false, rng)
    @test H isa SparseArrays.SparseMatrixCSC{Float64,Int64}
    @test norm(H - H') ≈ 0.0
    @test norm(H - transpose(H)) ≈ 0.0
    p = length(H.nzval) / N^2
    @test abs(p - d) < 0.1
    λ = eigvals(Array(H))
    @test λ isa Vector{Float64}
    sort!(λ)
    #@show relerr(λ[1], -ρ), relerr(λ[end], ρ)
    @test λ[1] ≊ -ρ
    @test λ[end] ≊ ρ

    ρ = 3.0
    d = 0.2
    H = random_matrix(N; density=d, hermitian=true, complex=false, spectral_radius=ρ, rng)
    @test H isa SparseArrays.SparseMatrixCSC{Float64,Int64}
    p = length(H.nzval) / N^2
    @test abs(p - d) < 0.1
    λ = sort(eigvals(Array(H)))
    #@show relerr(λ[1], -ρ), relerr(λ[end], ρ)
    @test λ[1] ≊ -ρ
    @test λ[end] ≊ ρ

    ρ = 3.0
    d = 0.2
    H = random_matrix(
        N;
        density=d,
        hermitian=true,
        complex=false,
        spectral_radius=ρ,
        rng,
        exact_spectral_radius=true
    )
    @test H isa SparseArrays.SparseMatrixCSC{Float64,Int64}
    λ = sort(eigvals(Array(H)))
    @test λ[1] ≈ -ρ
    @test λ[end] ≈ ρ

end


@testset "Random Hermitian Sparse Complex Matrix" begin

    H = random_matrix(N; density=0.5, hermitian=true)
    @test H isa SparseArrays.SparseMatrixCSC{ComplexF64,Int64}
    @test size(H) == (N, N)

    rng = StableRNG(3677556226)

    ρ = 1.0
    d = 0.5
    H = random_matrix(N; density=d, hermitian=true, rng)
    @test H isa SparseArrays.SparseMatrixCSC{ComplexF64,Int64}
    @test norm(H - H') ≈ 0.0
    @test norm(H - transpose(H)) > 0.0
    p = length(H.nzval) / N^2
    @test abs(p - d) < 0.1
    λ = eigvals(Array(H))
    @test λ isa Vector{Float64}
    sort!(λ)
    #@show relerr(λ[1], -ρ), relerr(λ[end], ρ)
    @test λ[1] ≊ -ρ
    @test λ[end] ≊ ρ

    ρ = 3.0
    d = 0.2
    H = random_matrix(N; density=d, hermitian=true, spectral_radius=ρ, rng)
    @test H isa SparseArrays.SparseMatrixCSC{ComplexF64,Int64}
    p = length(H.nzval) / N^2
    @test abs(p - d) < 0.1
    λ = sort(eigvals(Array(H)))
    #@show relerr(λ[1], -ρ), relerr(λ[end], ρ)
    @test λ[1] ≊ -ρ
    @test λ[end] ≊ ρ

    ρ = 3.0
    d = 0.2
    H = random_matrix(
        N;
        density=d,
        hermitian=true,
        spectral_radius=ρ,
        rng,
        exact_spectral_radius=true
    )
    @test H isa SparseArrays.SparseMatrixCSC{ComplexF64,Int64}
    λ = sort(eigvals(Array(H)))
    @test λ[1] ≈ -ρ
    @test λ[end] ≈ ρ

end


@testset "Random General Dense Real Matrix" begin

    H = random_matrix(N; complex=false)
    @test H isa Matrix{Float64}
    @test size(H) == (N, N)

    rng = StableRNG(1582783931)

    ρ = 1.0
    H = random_matrix(N; complex=false, rng)
    @test H isa Matrix{Float64}
    λ = eigvals(H)
    @test λ isa Vector{ComplexF64}
    @test maximum(imag.(λ)) > 0.5
    λ_max = maximum(abs.(λ))
    #@show relerr(λ_max, ρ)
    @test λ_max ≊ ρ

    ρ = 3.0
    H = random_matrix(N; spectral_radius=ρ, complex=false, rng, exact_spectral_radius=true)
    @test H isa Matrix{Float64}
    λ = eigvals(H)
    λ_max = maximum(abs.(λ))
    @test λ_max ≈ ρ

end


@testset "Random General Dense Complex Matrix" begin

    H = random_matrix(N)
    @test H isa Matrix{ComplexF64}
    @test maximum(imag.(H)) > 0.0
    @test size(H) == (N, N)

    rng = StableRNG(2137585613)

    ρ = 1.0
    H = random_matrix(N; rng)
    @test H isa Matrix{ComplexF64}
    λ = eigvals(H)
    @test λ isa Vector{ComplexF64}
    @test maximum(imag.(λ)) > 0.5
    λ_max = maximum(abs.(λ))
    #@show relerr(λ_max, ρ)
    @test λ_max ≊ ρ

    ρ = 3.0
    H = random_matrix(N; spectral_radius=ρ, rng, exact_spectral_radius=true)
    @test H isa Matrix{ComplexF64}
    λ = eigvals(H)
    λ_max = maximum(abs.(λ))
    @test λ_max ≈ ρ

end


@testset "Random General Sparse Real Matrix" begin

    H = random_matrix(N; density=0.5, complex=false)
    @test H isa SparseArrays.SparseMatrixCSC{Float64,Int64}
    @test size(H) == (N, N)

    rng = StableRNG(1710562897)

    ρ = 1.0
    d = 0.5
    H = random_matrix(N; density=d, complex=false, rng)
    @test H isa SparseArrays.SparseMatrixCSC{Float64,Int64}
    p = length(H.nzval) / N^2
    @test abs(p - d) < 0.1
    λ = eigvals(Array(H))
    @test λ isa Vector{ComplexF64}
    @test maximum(imag.(λ)) > 0.5
    λ_max = maximum(abs.(λ))
    #@show relerr(λ_max, ρ)
    @test λ_max ≊ ρ

    ρ = 3.0
    d = 0.2
    H = random_matrix(N; density=d, complex=false, spectral_radius=ρ, rng)
    @test H isa SparseArrays.SparseMatrixCSC{Float64,Int64}
    p = length(H.nzval) / N^2
    @test abs(p - d) < 0.1
    λ = eigvals(Array(H))
    @test λ isa Vector{ComplexF64}
    @test maximum(imag.(λ)) > 0.5
    λ_max = maximum(abs.(λ))
    #@show relerr(λ_max, ρ)
    @test λ_max ≊ ρ

    ρ = 3.0
    d = 0.2
    H = random_matrix(
        N;
        density=d,
        complex=false,
        spectral_radius=ρ,
        rng,
        exact_spectral_radius=true
    )
    @test H isa SparseArrays.SparseMatrixCSC{Float64,Int64}
    λ = eigvals(Array(H))
    λ_max = maximum(abs.(λ))
    @test λ_max ≈ ρ

end


@testset "Random General Sparse Complex Matrix" begin

    H = random_matrix(N; density=0.5)
    @test H isa SparseArrays.SparseMatrixCSC{ComplexF64,Int64}
    @test maximum(imag.(H)) > 0.0
    @test size(H) == (N, N)

    rng = StableRNG(2799774491)

    ρ = 1.0
    d = 0.5
    H = random_matrix(N; density=d, rng)
    @test H isa SparseArrays.SparseMatrixCSC{ComplexF64,Int64}
    p = length(H.nzval) / N^2
    @test abs(p - d) < 0.1
    λ = eigvals(Array(H))
    @test λ isa Vector{ComplexF64}
    @test maximum(imag.(λ)) > 0.5
    λ_max = maximum(abs.(λ))
    #@show relerr(λ_max, ρ)
    @test λ_max ≊ ρ

    ρ = 3.0
    d = 0.2
    H = random_matrix(N; density=d, spectral_radius=ρ, rng)
    @test H isa SparseArrays.SparseMatrixCSC{ComplexF64,Int64}
    p = length(H.nzval) / N^2
    @test abs(p - d) < 0.1
    λ = eigvals(Array(H))
    @test λ isa Vector{ComplexF64}
    @test maximum(imag.(λ)) > 0.5
    λ_max = maximum(abs.(λ))
    #@show relerr(λ_max, ρ)
    @test λ_max ≊ ρ

    ρ = 3.0
    d = 0.2
    H = random_matrix(N; density=d, spectral_radius=ρ, rng, exact_spectral_radius=true)
    @test H isa SparseArrays.SparseMatrixCSC{ComplexF64,Int64}
    λ = eigvals(Array(H))
    λ_max = maximum(abs.(λ))
    @test λ_max ≈ ρ

end


@testset "Random dynamic generator" begin


    tlist = collect(range(0, 100, length=1001))

    H = random_dynamic_generator(N, tlist)
    @test H isa Generator{Matrix{Float64},Vector{Float64}}
    @test length(H.ops) == 2

    state = random_state_vector(N)
    @test check_generator(H; state, tlist)

    rng = StableRNG(2316393754)

    H = random_dynamic_generator(N, tlist; rng)
    H_n = Array(evaluate(H, tlist, 1))
    @test size(H_n) == (N, N)
    λ = eigvals(H_n)
    @test λ isa Vector{Float64}
    @test -1 < λ[1] < 0
    @test 0 < λ[end] < 1

    H = random_dynamic_generator(N, tlist; rng, number_of_controls=3)
    @test H isa Generator{Matrix{Float64},Vector{Float64}}
    @test length(H.ops) == 4

    H = random_dynamic_generator(N, tlist; rng, complex=true)
    @test H isa Generator{Matrix{ComplexF64},Vector{Float64}}


    H = random_dynamic_generator(N, tlist; hermitian=false)
    @test H isa Generator{Matrix{Float64},Vector{Float64}}
    @test length(H.ops) == 2
    H_n = Array(evaluate(H, tlist, 1))
    λ = eigvals(H_n)
    @test λ isa Vector{ComplexF64}


    H = random_dynamic_generator(
        N,
        tlist;
        rng,
        hermitian=true,
        density=0.5,
        spectral_envelope=2.0
    )
    @test H isa Generator{SparseArrays.SparseMatrixCSC{Float64,Int64},Vector{Float64}}
    λ = reduce(vcat, [eigvals(Array(evaluate(H, tlist, n))) for n = 1:20:1000])
    @test λ isa Vector{Float64}
    @test -2 < λ[1] < -1
    @test 1 < λ[end] < 2

    H = random_dynamic_generator(
        N,
        tlist;
        rng,
        hermitian=true,
        density=0.5,
        spectral_envelope=2.0,
        exact_spectral_envelope=true
    )
    @test H isa Generator{SparseArrays.SparseMatrixCSC{Float64,Int64},Vector{Float64}}
    λ = vcat(
        eigvals(Array(evaluate(H, tlist, 1; vals_dict=IdDict(H.amplitudes[1] => 1.0)))),
        eigvals(Array(evaluate(H, tlist, 1; vals_dict=IdDict(H.amplitudes[1] => -1.0))))
    )
    @test λ isa Vector{Float64}
    @test abs(maximum(abs.(λ)) - 2.0) < 1e-5  # exact!
    λ = reduce(vcat, [eigvals(Array(evaluate(H, tlist, n))) for n = 1:50:1000])
    @test 1.9 < maximum(abs.(λ)) ≤ 2.0

    H = random_dynamic_generator(
        N,
        tlist;
        rng,
        hermitian=false,
        density=0.5,
        spectral_envelope=2.0,
        exact_spectral_envelope=true
    )
    @test H isa Generator{SparseArrays.SparseMatrixCSC{Float64,Int64},Vector{Float64}}
    λ = vcat(
        eigvals(Array(evaluate(H, tlist, 1; vals_dict=IdDict(H.amplitudes[1] => 1.0)))),
        eigvals(Array(evaluate(H, tlist, 1; vals_dict=IdDict(H.amplitudes[1] => -1.0))))
    )
    @test λ isa Vector{ComplexF64}
    @test abs(maximum(abs.(λ)) - 2.0) < 1e-5  # exact!
    λ = reduce(vcat, [eigvals(Array(evaluate(H, tlist, n))) for n = 1:50:1000])
    @test 1.9 < maximum(abs.(λ)) ≤ 2.0

end


@testset "Random state" begin

    Ψ = random_state_vector(N)
    @test Ψ isa Vector{ComplexF64}
    @test maximum(imag.(Ψ)) > 0.0
    @test norm(Ψ) ≈ 1.0

end
