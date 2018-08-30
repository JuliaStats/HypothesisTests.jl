using Documenter, HypothesisTests

makedocs(
    format = :html,
    sitename = "HypothesisTests.jl",
    modules = [HypothesisTests],
    pages = [
        "index.md",
        "methods.md",
        "parametric.md",
        "nonparametric.md",
        "time_series.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaStats/HypothesisTests.jl.git",
    target = "build",
    julia  = "1.0",
    deps = nothing,
    make = nothing
)
