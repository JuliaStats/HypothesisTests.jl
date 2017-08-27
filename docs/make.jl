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
    julia  = "0.6",
    deps = nothing,
    make = nothing
)
