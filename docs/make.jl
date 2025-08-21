using Documenter
using DocumenterCitations
using HypothesisTests

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

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
        "bibliography.md",
    ],
    checkdocs = :exports,
    plugins = [bib],
)

deploydocs(
    repo = "github.com/JuliaStats/HypothesisTests.jl.git",
    target = "build",
)
