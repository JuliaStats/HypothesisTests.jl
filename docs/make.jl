using Documenter, HypothesisTests

makedocs(
    modules = [HypothesisTests],
    sitename = "HypothesisTests.jl",
    pages = [
        "index.md",
        "methods.md",
        "parametric.md",
        "nonparametric.md",
        "time_series.md",
        "multivariate.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaStats/HypothesisTests.jl.git",
    target = "build",
)
