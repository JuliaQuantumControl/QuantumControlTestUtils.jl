module QuantumControlTestUtils

include("runner.jl")
export test
export show_coverage, generate_coverage_html

include("random.jl")
include("dummy_problem.jl")
include("logging.jl")

end
