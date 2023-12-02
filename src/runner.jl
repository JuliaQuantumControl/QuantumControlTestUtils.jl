# COV_EXCL_START
using Printf
using Logging
using Coverage
using PrettyTables: pretty_table, Highlighter
using LocalCoverage:
    LocalCoverage,
    eval_coverage_metrics,
    PackageCoverage,
    FileCoverageSummary,
    CoverageTools

function _coverage(metric::Union{PackageCoverage,FileCoverageSummary})
    return 100.0 * (float(metric.lines_hit) / metric.lines_tracked)
end

function format_line(metric::Union{PackageCoverage,FileCoverageSummary})
    name = (metric isa PackageCoverage ? "TOTAL" : relpath(metric.filename))
    lines_hit = @sprintf("%3d", metric.lines_hit)
    lines_tracked = @sprintf("%3d", metric.lines_tracked)
    lines_missing = @sprintf("%3d", metric.lines_tracked - metric.lines_hit)
    coverage = isnan(_coverage(metric)) ? "-" : @sprintf("%3.0f%%", _coverage(metric))
    return hcat(name, lines_tracked, lines_hit, lines_missing, coverage)
end


"""Collect coverage information from all `.cov` and `.info` files recursively
found in `root`. Return a Vector of `FileCoverage` objects filtered to files in
`paths` (relative to `root`).
"""
function collect_coverage(paths::Vector{String}=["src", "ext"]; root=pwd())
    root = abspath(root)
    local coverage
    logger = Logging.SimpleLogger(stderr, Logging.Error)
    Logging.with_logger(logger) do
        coverage = merge_coverage_counts(
            Coverage.process_folder(root),  # .cov files in root
            Coverage.LCOV.readfolder(root),  # tracefile.info
        )
    end
    coverage = filter(coverage) do covitem
        any(startswith(abspath(covitem.filename), joinpath(root, path)) for path in paths)
    end
    return coverage
end


"""Print out a coverage summary from existing coverage data.

```julia
show_coverage(paths=["./src", "./ext"]; root=pwd(), sort_by=nothing)
```

prints a a table showing the tracked files in `paths`, the total number of
tracked lines in that file ("Total"), the number of lines with coverage
("Hit"), the number of lines without coverage ("Missed") and the "Coverage" as
a percentage.

The coverage data is collected from `.cov` files in `paths` as well as
`tracefile-*.info` files in `root`.

Optionally, the table can be sorted by passing the name of a column to
`sort_by`, e..g. `sort_py=:Missed`.
"""
function show_coverage(paths=["src", "ext"]; root=pwd(), kwargs...)
    coverage = collect_coverage(paths; root=root)
    metrics = eval_coverage_metrics(coverage, root)
    show_coverage(metrics; kwargs...)
end

function show_coverage(metrics::PackageCoverage; sort_by=nothing)

    file_metrics = metrics.files
    sorter = Dict(
        :Total => (m -> m.lines_tracked),
        :Hit => (m -> m.lines_hit),
        :Missed => (m -> (m.lines_tracked - m.lines_hit)),
        :Coverage => (m -> _coverage(m)),
    )
    if !isnothing(sort_by)
        if sort_by ∈ keys(sorter)
            sort!(file_metrics; by=sorter[sort_by])
        else
            error("Cannot sort by $sort_by, must be one of $(keys(sorter))")
        end
    end

    table = reduce(vcat, map(format_line, [file_metrics..., metrics]))

    row_coverage = [[_coverage(m) for m in file_metrics]... _coverage(metrics)]

    highlighters = (
        Highlighter(
            (data, i, j) -> j == 5 && row_coverage[i] <= 50,
            bold=true,
            foreground=:red,
        ),
        Highlighter((data, i, j) -> j == 5 && row_coverage[i] >= 90, foreground=:green),
    )

    table_str = pretty_table(
        table,
        header=["File name", "Total", "Hit", "Missed", "Coverage"],
        alignment=[:l, :r, :r, :r, :r],
        crop=:none,
        linebreaks=true,
        columns_width=[maximum(length.(table[:, 1])), 6, 6, 6, 8],
        autowrap=false,
        highlighters=highlighters,
        body_hlines=[size(table, 1) - 1],
    )
    println(table_str)

end

_show_coverage_func = show_coverage


"""Generate an HTML report for existing coverage data.

```julia
generate_coverage_html(
    paths=["src", "ext"]; root=pwd(), covdir="coverage", genhtml="genhtml"
)
```

creates a folder `covdir` in `root` and use the external `genhtml` program to
write an HTML coverage report into that folder.
"""
function generate_coverage_html(paths=["src", "ext"]; root=pwd(), kwargs...)
    coverage = collect_coverage(paths; root=root)
    generate_coverage_html(root, coverage; kwargs...)
end


function generate_coverage_html(
    root::String,
    coverage::Vector{LocalCoverage.CoverageTools.FileCoverage};
    covdir="coverage",
    genhtml="genhtml"
)
    root = abspath(normpath(root))
    covdir = normpath(root, covdir)
    mkpath(covdir)
    tracefile = joinpath(covdir, "lcov.info")
    LocalCoverage.CoverageTools.LCOV.writefile(tracefile, coverage)
    branch = try
        strip(read(`git rev-parse --abbrev-ref HEAD`, String))
    catch e
        @warn "git branch could not be detected.\nError message: $(sprint(Base.showerror, e))"
    end
    title = isnothing(branch) ? "N/A" : "on branch $(branch)"
    try
        run(`$(genhtml) -t $(title) -o $(covdir) $(tracefile)`)
    catch e
        @error(
            "Failed to run $(genhtml). Check that lcov is installed.\nError message: $(sprint(Base.showerror, e))"
        )
    end
    @info(
        "Generated coverage HTML. Serve with 'LiveServer.serve(dir=\"$(relpath(covdir, pwd()))\")'"
    )
end


"""Run a package test-suite in a subprocess.

```julia
test(
    file="test/runtests.jl";
    root=pwd(),
    project="test",
    code_coverage=".coverage/tracefile-%p.info",
    show_coverage=(code_coverage != "none"),
    color=<inherit>,
    compiled_modules=<inherit>,
    startup_file=<inherit>,
    depwarn=<inherit>,
    inline=<inherit>,
    check_bounds="yes",
    track_allocation=<inherit>,
    threads=<inherit>,
    genhtml=false,
    covdir="coverage"
)
```

runs the test suite of the package located at `root` by running `include(file)`
inside a new julia process.

This is similar to what `Pkg.test()` does, but differs in the "sandboxing"
approach. While `Pkg.test()` creates a new temporary sandboxed environment,
`test()` uses an existing environment in `project` (the `test` subfolder by
default). This allows testing against the dev-versions of other packages. It
requires that the `test` folder contains both a `Project.toml` and a
`Manifest.toml` file.

The `test()` function also differs from directly including `test/runtests.jl`
in the REPL in that it can generate coverage data and reports (this is only
possible when running tests in a subprocess).

If `show_coverage` is passed as `true` (default), a coverage summary is shown.
Further, if `genhtml` is `true`, a full HTML coverage report will be generated
in `covdir` (relative to `root`). This requires the `genhtml` executable (part
of the [lcov](http://ltp.sourceforge.net/coverage/lcov.php) package). Instead
of `true`, it is also possible to pass the path to the `genhtml` executable.

All other keyword arguments correspond to the respective command line flag for
the `julia` executable that is run as the subprocess.

This function is intended to be exposed in a project's development-REPL.
"""
function test(
    file="test/runtests.jl";
    root=pwd(),
    project="test",
    code_coverage=joinpath(".coverage", "tracefile-%p.info"),
    show_coverage=(code_coverage != "none"),
    color=(Base.have_color === nothing ? "auto" : Base.have_color ? "yes" : "no"),
    compiled_modules=(Bool(Base.JLOptions().use_compiled_modules) ? "yes" : "no"),
    startup_file=(Base.JLOptions().startupfile == 1 ? "yes" : "no"),
    depwarn=(Base.JLOptions().depwarn == 2 ? "error" : "yes"),
    inline=(Bool(Base.JLOptions().can_inline) ? "yes" : "no"),
    track_allocation=(("none", "user", "all")[Base.JLOptions().malloc_log+1]),
    check_bounds="yes",
    threads=Threads.nthreads(),
    genhtml::Union{Bool,AbstractString}=false,
    covdir="coverage"
)
    root = abspath(normpath(root))
    if !(startswith(code_coverage, "@") || code_coverage in ["none", "user", "all"])
        code_coverage = normpath(root, code_coverage)
        if !endswith(code_coverage, ".info")
            @error "`code_coverage` must have `.info` extension" code_coverage
        end
        mkpath(dirname(code_coverage))
    end
    julia = Base.julia_cmd().exec[1]
    cmd = [
        julia,
        "--project=$project",
        "--color=$color",
        "--compiled-modules=$compiled_modules",
        "--startup-file=$startup_file",
        "--code-coverage=$code_coverage",
        "--track-allocation=$track_allocation",
        "--depwarn=$depwarn",
        "--check-bounds=$check_bounds",
        "--threads=$threads",
        "--inline=$inline",
        "--eval",
        "include(\"$file\")"
    ]
    @info "Running '$(join(cmd, " "))' in subprocess"
    run(Cmd(Cmd(cmd), dir=root))
    if show_coverage || genhtml
        logger = Logging.SimpleLogger(stderr, Logging.Error)
        coverage = collect_coverage(; root=root)
        if show_coverage
            metrics = eval_coverage_metrics(coverage, root)
            _show_coverage_func(metrics)
        end
        (genhtml === true) && (genhtml = "genhtml")
        (genhtml === false) && (genhtml = "")
        if !isempty(genhtml)
            generate_html_coverage(root; covdir, genhtml)
        end
    end
end
# COV_EXCL_STOP
