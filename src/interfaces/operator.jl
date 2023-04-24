using Test

using LinearAlgebra
using QuantumControl.Controls: get_controls, evaluate


"""Check that `op` is a valid operator that can be applied to `state`.

```julia
@test check_operator(op; state)
```

verifies the given `op` relative to `state`.

The specific requirements for `op` are:

* `op` must not be time-dependent: `evaluate(op, 0.0) == op`
* `op` must not contain any controls: `length(get_controls(op)) == 0`
* The 5-argument `mul!` must apply `op` to the given `state`
* `op * state` must be defined

"""
function check_operator(op; state)

    success = true

    Ψ = state

    @test check_state(state)
    ϕ0 = similar(state)
    ϕ = copy(ϕ0)

    try
        H = evaluate(op, 0.0)
        @test H == op
    catch exc
        @error "op must not be time-dependent: $exc"
        success = false
    end

    try
        @test length(get_controls(op)) == 0
    catch exc
        @error "op must not contain any controls: $exc"
        success = false
    end

    try
        @test mul!(ϕ, op, Ψ, 0.5, 0.5) ≡ ϕ
    catch exc
        @error "The 5-argument `mul!` must apply `op` to the given `state`: $exc"
        success = false
    end

    try
        ϕ1 = (op * Ψ)
        ϕ2 = 0.5 * ϕ1
        # check correctness of ϕ from mul!
        @test norm(ϕ - (0.5 * ϕ0 + ϕ2)) < 1e-12
    catch exc
        @error "`op * state` must be defined: $exc"
        success = false
    end

    # TODO: scalar multiplication?

    # TODO: dot?

    return success

end
