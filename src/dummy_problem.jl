using QuantumPropagators.Controls: get_controls, discretize, discretize_on_midpoints
using QuantumControl: Objective, ControlProblem


"""Set up a dummy control problem.

```julia
problem = dummy_control_problem(;
    N=10, n_objectives=1, n_controls=1, n_steps=50, dt=1.0, sparsity=0.5,
    complex_operators=true, hermitian=true, kwargs...)
```

Sets up a control problem with random (sparse) Hermitian matrices.

# Arguments

* `N`: The dimension of the Hilbert space
* `n_objectives`: The number of objectives in the optimization. All objectives
  will have the same Hamiltonian, but random initial and target states.
* `n_controls`: The number of controls, that is, the number of control terms in
  the control Hamiltonian. Each control is an array of random values,
  normalized on the intervals of the time grid.
* `n_steps`: The number of time steps (intervals of the time grid)
* `dt`: The time step
* `sparsity`: The sparsity of the Hamiltonians, as a number between 0.0 and
  1.0. For `sparsity=1.0`, the Hamiltonians will be dense matrices.
* `complex_operators`: Whether or not the drift/control operators will be
  complex-valued or real-valued.
* `hermitian`: Whether or not all drift/control operators will be Hermitian matrices.
* `kwargs`: All other keyword arguments are passed on to
  [`ControlProblem`](@ref)
"""
function dummy_control_problem(;
    N=10,
    n_objectives=1,
    n_controls=1,
    n_steps=50,
    dt=1.0,
    sparsity=0.5,
    complex_operators=true,
    hermitian=true,
    kwargs...
)

    tlist = collect(range(0; length=(n_steps + 1), step=dt))
    pulses = [rand(length(tlist) - 1) for l = 1:n_controls]
    for l = 1:n_controls
        # we normalize on the *intervals*, not on the time grid points
        pulses[l] ./= norm(pulses[l])
    end
    controls = [discretize(pulse, tlist) for pulse in pulses]

    function random_op(N, ρ, sparsity, complex_operators, hermitian)
        if sparsity < 1.0
            if hermitian
                if complex_operators
                    H = random_hermitian_sparse_matrix(N, ρ, sparsity)
                end
            else
                if complex_operators
                    H = random_complex_sparse_matrix(N, ρ, sparsity)
                end
            end
        else
            if hermitian
                if complex_operators
                    H = random_hermitian_matrix(N, ρ)
                else
                    H = random_hermitian_real_matrix(N, ρ)
                end
            else
                if complex_operators
                    H = random_complex_matrix(N, ρ)
                else
                    H = random_real_matrix(N, ρ)
                end
            end
        end
    end

    hamiltonian = []
    H_0 = random_op(N, 1.0, sparsity, complex_operators, hermitian)
    push!(hamiltonian, H_0)
    for control ∈ controls
        H_c = random_op(N, 1.0, sparsity, complex_operators, hermitian)
        push!(hamiltonian, (H_c, control))
    end

    objectives = [
        Objective(;
            initial_state=random_state_vector(N),
            generator=tuple(hamiltonian...),
            target_state=random_state_vector(N)
        ) for k = 1:n_objectives
    ]

    return ControlProblem(
        objectives=objectives,
        pulse_options=Dict(
            control => Dict(:lambda_a => 1.0, :update_shape => t -> 1.0) for
            control in controls
        ),
        tlist=tlist,
        kwargs...
    )
end


"""Result returned by [`optimize_with_dummy_method`](@ref)."""
mutable struct DummyOptimizationResult
    tlist::Vector{Float64}
    iter_start::Int64  # the starting iteration number
    iter_stop::Int64 # the maximum iteration number
    iter::Int64  # the current iteration number
    J_T::Float64  # the current value of the final-time functional J_T
    J_T_prev::Float64  # previous value of J_T
    guess_controls::Vector{Vector{Float64}}
    optimized_controls::Vector{Vector{Float64}}
    converged::Bool
    message::String

    function DummyOptimizationResult(problem)
        tlist = problem.tlist
        controls = get_controls(problem.objectives)
        iter_start = get(problem.kwargs, :iter_start, 0)
        iter = iter_start
        iter_stop = get(problem.kwargs, :iter_stop, 20)
        guess_controls = [discretize(control, tlist) for control in controls]
        J_T = 0.0
        J_T_prev = 0.0
        optimized_controls = [copy(guess) for guess in guess_controls]
        converged = false
        message = "in progress"
        new(
            tlist,
            iter_start,
            iter_stop,
            iter,
            J_T,
            J_T_prev,
            guess_controls,
            optimized_controls,
            converged,
            message
        )
    end

end

struct DummyOptimizationWrk
    objectives
    adjoint_objectives
    kwargs
    controls
    pulses0::Vector{Vector{Float64}}
    pulses1::Vector{Vector{Float64}}
    result
end


function DummyOptimizationWrk(problem)
    objectives = [obj for obj in problem.objectives]
    adjoint_objectives = [adjoint(obj) for obj in problem.objectives]
    controls = get_controls(objectives)
    kwargs = Dict(problem.kwargs)
    tlist = problem.tlist
    if haskey(kwargs, :continue_from)
        @info "Continuing previous optimization"
        result = kwargs[:continue_from]
        if !(result isa DummyOptimizationResult)
            result = convert(DummyOptimizationResult, result)
        end
        result.iter_stop = get(problem.kwargs, :iter_stop, 20)
        result.converged = false
        result.message = "in progress"
        pulses0 = [
            discretize_on_midpoints(control, tlist) for control in result.optimized_controls
        ]
    else
        result = DummyOptimizationResult(problem)
        pulses0 = [discretize_on_midpoints(control, tlist) for control in controls]
    end
    pulses1 = [copy(pulse) for pulse in pulses0]
    return DummyOptimizationWrk(
        objectives,
        adjoint_objectives,
        kwargs,
        controls,
        pulses0,
        pulses1,
        result,
    )
end

function update_result!(wrk::DummyOptimizationWrk, ϵ⁽ⁱ⁺¹⁾, i::Int64)
    res = wrk.result
    res.J_T_prev = res.J_T
    res.J_T = sum([norm(ϵ) for ϵ ∈ ϵ⁽ⁱ⁺¹⁾])
    (i > 0) && (res.iter = i)
    if i >= res.iter_stop
        res.converged = true
        res.message = "Reached maximum number of iterations"
    end
end


"""Run a dummy optimization.

```julia
result = optimize(problem, method=:dummymethod)
```

runs through and "optimization" of the given `problem` where in each iteration,
the amplitude of the guess pulses is diminished by 10%. The (summed) vector
norm of the the control serves as the value of the optimization functional.
"""
function optimize_with_dummy_method(problem)
    # This is connected to the main `optimize` method in the main
    # QuantumcontrolBase.jl
    iter_start = get(problem.kwargs, :iter_start, 0)
    check_convergence! = get(problem.kwargs, :check_convergence, res -> res)
    wrk = DummyOptimizationWrk(problem)
    ϵ⁽ⁱ⁾ = wrk.pulses0
    ϵ⁽ⁱ⁺¹⁾ = wrk.pulses1
    update_result!(wrk, ϵ⁽ⁱ⁺¹⁾, 0)
    println("# iter\tJ_T")
    @printf("%6d\t%.2e\n", 0, wrk.result.J_T)
    i = wrk.result.iter  # = 0, unless continuing from previous optimization
    while !wrk.result.converged
        i = i + 1
        for l = 1:length(ϵ⁽ⁱ⁺¹⁾)
            ϵ⁽ⁱ⁺¹⁾[l] .= 0.9 * ϵ⁽ⁱ⁾[l]
        end
        update_result!(wrk, ϵ⁽ⁱ⁺¹⁾, i)
        @printf("%6d\t%.2e\n", i, wrk.result.J_T)
        check_convergence!(wrk.result)
        ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾ = ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾
    end
    for l = 1:length(ϵ⁽ⁱ⁾)
        wrk.result.optimized_controls[l] = discretize(ϵ⁽ⁱ⁾[l], problem.tlist)
    end
    return wrk.result
end
