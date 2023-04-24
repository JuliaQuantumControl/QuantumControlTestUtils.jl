using Test

using LinearAlgebra


"""Check that `state` is valid.

```julia
@test check_state(state)
```

verifies the following requirements:

* `similar(state)` must be defined
* The inner product of two states must be a complex number
* The `norm` of `state` must be defined via the inner product.
* States must be able to be added and subtracted
* `copy(state)` must be defined
* `c * state` for a scalar `c` must be defined
"""
function check_state(state)

    success = true

    try
        @test similar(state) ≢ nothing
    catch exc
        @error "similar(state) must be defined: $exc"
        success = false
    end

    try
        @test (state ⋅ state) isa Complex
    catch exc
        @error "the inner product of two states must be a complex number: $exc"
        success = false
    end

    try
        @test norm(state) ≈ sqrt(state ⋅ state)
    catch exc
        @error "The norm of a state must be defined via the inner product: $exc"
        success = false
    end

    try
        @test norm(state - state) ≈ 0.0
        @test norm(state + state) ≤ 2 * norm(state)
    catch exc
        @error "States must be able to be added and subtracted: $exc"
        success = false
    end

    try
        ϕ = copy(state)
        @test norm(ϕ - state) ≈ 0.0
    catch exc
        @error "copy(state) must be defined: $exc"
        success = false
    end

    try
        ϕ = 0.5 * state
        @test norm(ϕ) < norm(state)
    catch exc
        @error "`c * state` for a scalar `c` must be defined: $exc"
        success = false
    end

    # TODO: for_inplace
    # * copyto!
    # * lmul!
    # * axpy!

    return success

end
