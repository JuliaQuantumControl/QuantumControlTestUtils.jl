module QuantumControlTestUtils

include("runner.jl")
export test
export show_coverage, generate_coverage_html

include("random.jl")  # submodule RandomObjects

include("dummy_problem.jl")  # submodule DummyOptimization

include("logging.jl")

end
