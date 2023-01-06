module RandomObjects

using Random
using LinearAlgebra
using SparseArrays

export random_state_vector, random_matrix


"""Construct a random matrix.

```julia
Ĥ = random_matrix(N; kwargs...)
```

by default initializes `Ĥ` as a general complex ``N×N`` matrix with a spectral
radius of approximately 1.0. Keyword arguments allow to initialize real or
complex, Hermitian or non-Hermitian, dense or sparse matrices with arbitrary
spectral radius. The non-zero entries in Ĥ will be uniformly distributed around
zero, with a range of values that depends on `N` and the desired spectral
radius.

# Keyword arguments

* `density=1.0`: A number > 0.0 and ≤ 1.0. Any value < 1.0 implies a sparse
  matrix where `density` is the approximate fraction of non-zero elements to
  total elements
* `complex=true`: Whether the matrix should be complex-valued (default) or
  real-valued
* `hermitian=false`: Whether the matrix should be general (default) or
  Hermitian (real eigenvalues)
* `spectral_radius=1.0`: The approximate spectral radius, i.e. maximum absolute
  eigenvalue. This is according to [Girko-Ginibri's circular
  law](https://www.johndcook.com/blog/2018/07/27/circular-law/), in the limit
  of large ``N``
* `exact_spectral_radius=false`: If given as `true`, ensure that the
  `spectral_radius` is exact. This is done via diagonalization, so it is only
  feasible for moderately large dimensions N. On the other hand, for large `N`,
  the `spectral_radius`, respectively the circular law becomes more exact
  anyway.
* `rng=Random.GLOBAL_RNG`: The random number generator to use. The call
  `Random.rand(rng, N, N)` must produces a real-valued ``N×N`` matrix with
  elements uniformly distributed between 0 and 1
"""
function random_matrix(
    N;
    density=1.0,
    complex=true,
    hermitian=false,
    spectral_radius=1.0,
    exact_spectral_radius=false,
    rng=Random.GLOBAL_RNG
)
    if (density ≤ 0.0) || (density > 1.0)
        error("density must be in (0, 1]")
    end
    if complex
        if density < 1.0  # sparse matrix
            if hermitian
                random_hermitian_sparse_matrix(
                    N,
                    spectral_radius,
                    density;
                    rng,
                    exact_spectral_radius
                )
            else
                random_complex_sparse_matrix(
                    N,
                    spectral_radius,
                    density;
                    rng,
                    exact_spectral_radius
                )
            end
        else  # dense matrix
            if hermitian
                random_hermitian_matrix(N, spectral_radius; rng, exact_spectral_radius)
            else
                random_complex_matrix(N, spectral_radius; rng, exact_spectral_radius)
            end
        end
    else  # real-valued matrix
        if density < 1.0  # sparse matrix
            if hermitian
                random_hermitian_sparse_real_matrix(
                    N,
                    spectral_radius,
                    density;
                    rng,
                    exact_spectral_radius
                )
            else
                random_real_sparse_matrix(
                    N,
                    spectral_radius,
                    density;
                    rng,
                    exact_spectral_radius
                )
            end
        else  # dense matrix
            if hermitian
                random_hermitian_real_matrix(N, spectral_radius; rng, exact_spectral_radius)
            else
                random_real_matrix(N, spectral_radius; rng, exact_spectral_radius)
            end
        end
    end
end


"""Construct a random complex matrix of size N×N with spectral radius ρ.

```julia
random_complex_matrix(N, ρ)
```
"""
function random_complex_matrix(N, ρ; rng=Random.GLOBAL_RNG, exact_spectral_radius=false)
    Δ = √(12 / N)
    X = Δ * (rand(rng, N, N) .- 0.5)
    Y = Δ * (rand(rng, N, N) .- 0.5)
    H = ρ * (X + Y * 1im) / √2
    if exact_spectral_radius
        λ = eigvals(H)
        Δ = maximum(abs.(λ))
        return (ρ / Δ) * H
    else
        return H
    end
end


"""Construct a random real-valued matrix of size N×N with spectral radius ρ.

```julia
random_real_matrix(N, ρ)
```
"""
function random_real_matrix(N, ρ; rng=Random.GLOBAL_RNG, exact_spectral_radius=false)
    Δ = √(12 / N)
    X = Δ * (rand(rng, N, N) .- 0.5)
    H = ρ * X
    if exact_spectral_radius
        λ = eigvals(H)
        Δ = maximum(abs.(λ))
        return (ρ / Δ) * H
    else
        return H
    end
end


"""Construct a random Hermitian complex matrix of size N×N with spectral radius ρ.

```julia
random_hermitian_matrix(N, ρ)
```
"""
function random_hermitian_matrix(N, ρ; rng=Random.GLOBAL_RNG, exact_spectral_radius=false)
    Δ = √(12 / N)
    X = Δ * (rand(rng, N, N) .- 0.5)
    Y = Δ * (rand(rng, N, N) .- 0.5)
    Z = (X + Y * 1im) / √2
    H = ρ * (Z + Z') / (2 * √2)
    if exact_spectral_radius
        λ = eigvals(H)
        λ₀ = λ[1]
        Δ = λ[end] - λ₀
        H_norm = (H - λ₀ * I) / Δ
        return 2 * ρ * H_norm - ρ * I
    else
        return H
    end
end


"""Construct a random Hermitian real matrix of size N×N with spectral radius ρ.

```julia
random_hermitian_real_matrix(N, ρ)
```
"""
function random_hermitian_real_matrix(
    N,
    ρ;
    rng=Random.GLOBAL_RNG,
    exact_spectral_radius=false
)
    Δ = √(12 / N)
    X = Δ * (rand(N, N) .- 0.5)
    H = ρ * (X + X') / (2 * √2)
    if exact_spectral_radius
        λ = eigvals(H)
        λ₀ = λ[1]
        Δ = λ[end] - λ₀
        H_norm = (H - λ₀ * I) / Δ
        return 2 * ρ * H_norm - ρ * I
    else
        return H
    end
end


"""Construct a random sparse complex matrix.

```julia
random_complex_sparse_matrix(N, ρ, density)
```

returns a matrix of size N×N with spectral radius ρ and the given density
(number between zero and one that is the approximate fraction of non-zero
elements).
"""
function random_complex_sparse_matrix(
    N,
    ρ,
    density;
    rng=Random.GLOBAL_RNG,
    exact_spectral_radius=false
)
    p = 1 - √(1 - density)
    Δ = √(12 / (p * N))
    X = sprand(rng, N, N, p)
    X.nzval .= Δ .* (X.nzval .- 0.5)
    Y = sprand(rng, N, N, p)
    Y.nzval .= Δ .* (Y.nzval .- 0.5)
    H = ρ * (X + Y * 1im) / √2
    if exact_spectral_radius
        λ = eigvals(Array(H))
        Δ = maximum(abs.(λ))
        return (ρ / Δ) * H
    else
        return H
    end
end


"""Construct a random sparse real-valued matrix.

```julia
random_real_sparse_matrix(N, ρ, density)
```

returns a matrix of size N×N with spectral radius ρ and the given density
(number between zero and one that is the approximate fraction of non-zero
elements).
"""
function random_real_sparse_matrix(
    N,
    ρ,
    density;
    rng=Random.GLOBAL_RNG,
    exact_spectral_radius=false
)
    p = density
    Δ = √(12 / (p * N))
    X = sprand(rng, N, N, density)
    X.nzval .= Δ .* (X.nzval .- 0.5)
    H = ρ * X
    if exact_spectral_radius
        λ = eigvals(Array(H))
        Δ = maximum(abs.(λ))
        return (ρ / Δ) * H
    else
        return H
    end
end


"""Construct a random sparse Hermitian matrix.

```julia
random_hermitian_sparse_matrix(N, ρ, density)
```

returns a matrix of size N×N with spectral radius ρ and the given density
(number between zero and one that is the approximate fraction of non-zero
elements).
"""
function random_hermitian_sparse_matrix(
    N,
    ρ,
    density;
    rng=Random.GLOBAL_RNG,
    exact_spectral_radius=false
)
    p = 1 - √(1 - density)
    Δ = √(12 / (p * N))
    X = sprand(rng, N, N, p)
    X.nzval .= Δ .* (X.nzval .- 0.5)
    Y = copy(X)
    Y.nzval .= Δ * (rand(rng, length(Y.nzval)) .- 0.5)
    Z = (X + Y * 1im) / √2
    H = ρ * (Z + Z') / (2 * √2)
    if exact_spectral_radius
        λ = eigvals(Array(H))
        λ₀ = λ[1]
        Δ = λ[end] - λ₀
        H_norm = (H - λ₀ * I) / Δ
        return 2 * ρ * H_norm - ρ * I
    else
        return H
    end
end


"""Construct a random sparse Hermitian real matrix.

```julia
random_hermitian_sparse_real_matrix(N, ρ, density)
```

returns a matrix of size N×N with spectral radius ρ and the given density
(number between zero and one that is the approximate fraction of non-zero
elements).
"""
function random_hermitian_sparse_real_matrix(
    N,
    ρ,
    density;
    rng=Random.GLOBAL_RNG,
    exact_spectral_radius=false
)
    p = 1 - √(1 - density)
    Δ = √(12 / (p * N))
    X = sprand(rng, N, N, p)
    X.nzval .= Δ .* (X.nzval .- 0.5)
    H = ρ * (X + X') / (2 * √2)
    if exact_spectral_radius
        λ = eigvals(Array(H))
        λ₀ = λ[1]
        Δ = λ[end] - λ₀
        H_norm = (H - λ₀ * I) / Δ
        return 2 * ρ * H_norm - ρ * I
    else
        return H
    end
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
