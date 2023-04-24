module QuantumControlTestUtils

include("runner.jl")
export test
export show_coverage, generate_coverage_html

include("random.jl")  # submodule RandomObjects

include("dummy_problem.jl")  # submodule DummyOptimization


module Interfaces

export check_operator, check_state
include(joinpath("interfaces", "state.jl"))
include(joinpath("interfaces", "operator.jl"))
# TODO: amplitude interface
# TODO: control interface
# TODO: generator interface
# TODO: propagator interface


end

include("logging.jl")

end
