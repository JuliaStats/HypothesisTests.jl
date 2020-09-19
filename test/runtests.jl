using Test, Random, Statistics

# Confidence intervals should be tested component-wise against
# expected results; splatting them into arrays and checking for
# approximate array equality is incorrect when either limit is
# infinite.
macro test_ci_approx(x::Expr, y::Expr)
    return quote
        Test.@test typeof($(esc(x))) <: Tuple{Real,Real}
        Test.@test typeof($(esc(x))) == typeof($(esc(y)))
        Test.@test all(map(isapprox, $(esc(x)), $(esc(y))))
    end
end

include("anderson_darling.jl")
include("augmented_dickey_fuller.jl")
include("bartlett.jl")
include("binomial.jl")
include("box_test.jl")
include("breusch_godfrey.jl")
include("circular.jl")
include("common.jl")
include("durbin_watson.jl")
include("fisher.jl")
include("jarque_bera.jl")
include("hotelling.jl")
include("kolmogorov_smirnov.jl")
include("kruskal_wallis.jl")
include("mann_whitney.jl")
include("correlation.jl")
include("permutation.jl")
include("power_divergence.jl")
include("show.jl")
include("t.jl")
include("wald_wolfowitz.jl")
include("wilcoxon.jl")
include("z.jl")
include("f.jl")
include("diebold_mariano.jl")
include("clark_west.jl")
include("white.jl")
include("var_equality.jl")
