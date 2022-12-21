module RandomObjects

using Random
using LinearAlgebra
using SparseArrays

export random_state_vector, random_matrix


function random_matrix(
    N;
    density=1.0,
    complex=true,
    hermitian=false,
    spectral_radius=1.0,
    rng=Random.GLOBAL_RNG
)
    if (density ≤ 0.0) || (density > 1.0)
        error("density must be in (0, 1]")
    end
    if complex
        if density < 1.0  # sparse matrix
            if hermitian
                random_hermitian_sparse_matrix(N, spectral_radius, density; rng)
            else
                random_complex_sparse_matrix(N, spectral_radius, density; rng)
            end
        else  # dense matrix
            if hermitian
                random_hermitian_matrix(N, spectral_radius; rng)
            else
                random_complex_matrix(N, spectral_radius; rng)
            end
        end
    else  # real-valued matrix
        if density < 1.0  # sparse matrix
            if hermitian
                random_hermitian_sparse_real_matrix(N, spectral_radius, density; rng)
            else
                random_real_sparse_matrix(N, spectral_radius, density; rng)
            end
        else  # dense matrix
            if hermitian
                random_hermitian_real_matrix(N, spectral_radius; rng)
            else
                random_real_matrix(N, spectral_radius; rng)
            end
        end
    end
end


"""Construct a random complex matrix of size N×N with spectral radius ρ.

```julia
random_complex_matrix(N, ρ)
```
"""
function random_complex_matrix(N, ρ; rng=Random.GLOBAL_RNG)
    Δ = √(12 / N)
    X = Δ * (rand(rng, N, N) .- 0.5)
    Y = Δ * (rand(rng, N, N) .- 0.5)
    H = ρ * (X + Y * 1im) / √2
end


"""Construct a random real-valued matrix of size N×N with spectral radius ρ.

```julia
random_real_matrix(N, ρ)
```
"""
function random_real_matrix(N, ρ; rng=Random.GLOBAL_RNG)
    Δ = √(12 / N)
    X = Δ * (rand(rng, N, N) .- 0.5)
    H = ρ * X
end


"""Construct a random Hermitian matrix of size N×N with spectral radius ρ.

```julia
random_hermitian_matrix(N, ρ)
```
"""
function random_hermitian_matrix(N, ρ; rng=Random.GLOBAL_RNG)
    Δ = √(12 / N)
    X = Δ * (rand(rng, N, N) .- 0.5)
    Y = Δ * (rand(rng, N, N) .- 0.5)
    Z = (X + Y * 1im) / √2
    H = ρ * (Z + Z') / (2 * √2)
end


"""Construct a random Hermitian real matrix of size N×N with spectral radius ρ.

```julia
random_hermitian_real_matrix(N, ρ)
```
"""
function random_hermitian_real_matrix(N, ρ; rng=Random.GLOBAL_RNG)
    Δ = √(12 / N)
    X = Δ * (rand(N, N) .- 0.5)
    H = ρ * (X + X') / (2 * √2)
end


"""Construct a random sparse complex matrix.

```julia
random_complex_sparse_matrix(N, ρ, density)
```

returns a matrix of size N×N with spectral radius ρ and the given density
(number between zero and one that is the approximate fraction of non-zero
elements).
"""
function random_complex_sparse_matrix(N, ρ, density; rng=Random.GLOBAL_RNG)
    p = 1 - √(1 - density)
    Δ = √(12 / (p * N))
    X = sprand(rng, N, N, p)
    X.nzval .= Δ .* (X.nzval .- 0.5)
    Y = sprand(rng, N, N, p)
    Y.nzval .= Δ .* (Y.nzval .- 0.5)
    H = ρ * (X + Y * 1im) / √2
end


"""Construct a random sparse real-valued matrix.

```julia
random_real_sparse_matrix(N, ρ, density)
```

returns a matrix of size N×N with spectral radius ρ and the given density
(number between zero and one that is the approximate fraction of non-zero
elements).
"""
function random_real_sparse_matrix(N, ρ, density; rng=Random.GLOBAL_RNG)
    p = density
    Δ = √(12 / (p * N))
    X = sprand(rng, N, N, density)
    X.nzval .= Δ .* (X.nzval .- 0.5)
    H = ρ * X
end


"""Construct a random sparse Hermitian matrix.

```julia
random_hermitian_sparse_matrix(N, ρ, density)
```

returns a matrix of size N×N with spectral radius ρ and the given density
(number between zero and one that is the approximate fraction of non-zero
elements).
"""
function random_hermitian_sparse_matrix(N, ρ, density; rng=Random.GLOBAL_RNG)
    p = 1 - √(1 - density)
    Δ = √(12 / (p * N))
    X = sprand(rng, N, N, p)
    X.nzval .= Δ .* (X.nzval .- 0.5)
    Y = copy(X)
    Y.nzval .= Δ * (rand(rng, length(Y.nzval)) .- 0.5)
    Z = (X + Y * 1im) / √2
    return ρ * (Z + Z') / (2 * √2)
end


"""Construct a random sparse Hermitian real matrix.

```julia
random_hermitian_sparse_real_matrix(N, ρ, density)
```

returns a matrix of size N×N with spectral radius ρ and the given density
(number between zero and one that is the approximate fraction of non-zero
elements).
"""
function random_hermitian_sparse_real_matrix(N, ρ, density; rng=Random.GLOBAL_RNG)
    p = 1 - √(1 - density)
    Δ = √(12 / (p * N))
    X = sprand(rng, N, N, p)
    X.nzval .= Δ .* (X.nzval .- 0.5)
    return ρ * (X + X') / (2 * √2)
end


"""Return a random, normalized Hilbert space state vector of dimension `N`.

```julia
random_state_vector(N; rng=GLOBAL_RNG)
```
"""
function random_state_vector(N; rng=Random.GLOBAL_RNG)
    Ψ = rand(rng, N) .* exp.((2π * im) .* rand(rng, N))
    Ψ ./= norm(Ψ)
    return Ψ
end


# undocumented, but useful for interactively getting RNG seeds in tests
randseed() = Int(UInt32(Random.bitrand(32).chunks[1]))

end
