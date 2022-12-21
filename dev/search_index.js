var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = QuantumControlTestUtils","category":"page"},{"location":"#QuantumControlTestUtils","page":"Home","title":"QuantumControlTestUtils","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for QuantumControlTestUtils.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [QuantumControlTestUtils]","category":"page"},{"location":"#QuantumControlTestUtils.DummyOptimizationResult","page":"Home","title":"QuantumControlTestUtils.DummyOptimizationResult","text":"Result returned by optimize_with_dummy_method.\n\n\n\n\n\n","category":"type"},{"location":"#QuantumControlTestUtils.QuantumTestLogger","page":"Home","title":"QuantumControlTestUtils.QuantumTestLogger","text":"Logger that stores all messages in a messages field (array of strings).\n\nThe in operator can be used to check whether a string or regex is in any of the recorded messages.\n\ntest_logger = QuantumTestLogger()\nwith_logger(test_logger) do\n    #...\nend\n@test \"message\" in logger\n\nThis is a simplified version of Logging.TestLogger.\n\n\n\n\n\n","category":"type"},{"location":"#QuantumControlTestUtils.dummy_control_problem-Tuple{}","page":"Home","title":"QuantumControlTestUtils.dummy_control_problem","text":"Set up a dummy control problem.\n\nproblem = dummy_control_problem(;\n    N=10, n_objectives=1, n_controls=1, n_steps=50, dt=1.0, sparsity=0.5,\n    complex_operators=true, hermitian=true, kwargs...)\n\nSets up a control problem with random (sparse) Hermitian matrices.\n\nArguments\n\nN: The dimension of the Hilbert space\nn_objectives: The number of objectives in the optimization. All objectives will have the same Hamiltonian, but random initial and target states.\nn_controls: The number of controls, that is, the number of control terms in the control Hamiltonian. Each control is an array of random values, normalized on the intervals of the time grid.\nn_steps: The number of time steps (intervals of the time grid)\ndt: The time step\nsparsity: The sparsity of the Hamiltonians, as a number between 0.0 and 1.0. For sparsity=1.0, the Hamiltonians will be dense matrices.\ncomplex_operators: Whether or not the drift/control operators will be complex-valued or real-valued.\nhermitian: Whether or not all drift/control operators will be Hermitian matrices.\nkwargs: All other keyword arguments are passed on to ControlProblem\n\n\n\n\n\n","category":"method"},{"location":"#QuantumControlTestUtils.generate_coverage_html","page":"Home","title":"QuantumControlTestUtils.generate_coverage_html","text":"Generate an HTML report for existing coverage data.\n\ngenerate_coverage_html(path=\"./\"; covdir=\"coverage\", genhtml=\"genhtml\")\n\ncreates a folder covdir and use the external genhtml program to write an HTML coverage report into that folder.\n\n\n\n\n\n","category":"function"},{"location":"#QuantumControlTestUtils.optimize_with_dummy_method-Tuple{Any}","page":"Home","title":"QuantumControlTestUtils.optimize_with_dummy_method","text":"Run a dummy optimization.\n\nresult = optimize(problem, method=:dummymethod)\n\nruns through and \"optimization\" of the given problem where in each iteration, the amplitude of the guess pulses is diminished by 10%. The (summed) vector norm of the the control serves as the value of the optimization functional.\n\n\n\n\n\n","category":"method"},{"location":"#QuantumControlTestUtils.show_coverage","page":"Home","title":"QuantumControlTestUtils.show_coverage","text":"Print out a coverage summary from existing coverage data.\n\nshow_coverage(path=\"./src\"; sort_by=nothing)\n\nprints a a table showing the tracked files in path, the total number of tracked lines in that file (\"Total\"), the number of lines with coverage (\"Hit\"), the number of lines without coverage (\"Missed\") and the \"Coverage\" as a percentage.\n\nThe coverage data is collected from .cov files in path.\n\nOptionally, the table can be sorted by passing the name of a column to sort_by, e..g. sort_py=:Missed.\n\n\n\n\n\n","category":"function"},{"location":"#QuantumControlTestUtils.test","page":"Home","title":"QuantumControlTestUtils.test","text":"Run a package test-suite in a subprocess.\n\ntest(\n    file=\"test/runtests.jl\";\n    root=pwd(),\n    project=\"test\",\n    code_coverage=\"user\",\n    show_coverage=(code_coverage == \"user\"),\n    color=<inherit>,\n    compiled_modules=<inherit>,\n    startup_file=<inherit>,\n    depwarn=<inherit>,\n    inline=<inherit>,\n    check_bounds=\"yes\",\n    track_allocation=<inherit>,\n    threads=<inherit>,\n    genhtml=false,\n    covdir=\"coverage\"\n)\n\nruns the test suite of the package located at root by running include(file) inside a new julia process.\n\nThis is similar to what Pkg.test() does, but differs in the \"sandboxing\" approach. While Pkg.test() creates a new temporary sandboxed environment, test() uses an existing environment in project (the test subfolder by default). This allows testing against the dev-versions of other packages. It requires that the test folder contains both a Project.toml and a Manifest.toml file.\n\nThe test() function also differs from directly including test/runtests.jl in the REPL in that it can generate coverage data and reports (this is only possible when running tests in a subprocess).\n\nIf show_coverage is passed as true (default), a coverage summary is shown. Further, if genhtml is true, a full HTML coverage report will be generated in covdir (relative to root). This requires the genhtml executable (part of the lcov package). Instead of true, it is also possible to pass the path to the genhtml exectuable.\n\nAll other keyword arguments correspond to the respective command line flag for the julia executable that is run as the subprocess.\n\nThis function is intended to be exposed in a project's development-REPL.\n\n\n\n\n\n","category":"function"}]
}
