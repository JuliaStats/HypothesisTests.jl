using Documenter
using DocumenterCitations
using HypothesisTests

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(
    modules = [HypothesisTests],
    sitename = "HypothesisTests.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        # The specification pages in the appendix run past the default 100 KiB advisory
        # limit; that is the shape of the page rather than a problem with it. The hard
        # threshold at 200 KiB is left where it is.
        size_threshold_warn = 150 * 1024,
    ),
    pages = [
        "index.md",
        "methods.md",
        "parametric.md",
        "nonparametric.md",
        "time_series.md",
        "multivariate.md",
        "Appendix: mathematical specifications" => [
            "specs/index.md",
            "specs/t_test.md",
            "specs/rank_tests.md",
        ],
        "bibliography.md",
    ],
    checkdocs = :exports,
    plugins = [bib],
)

deploydocs(
    repo = "github.com/JuliaStats/HypothesisTests.jl.git",
    target = "build",
)
