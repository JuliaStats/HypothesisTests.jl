# D'Agostino test of skewness

export AgostinoTest

struct AgostinoTest <: HypothesisTest
    nobs::Int               # number of observations
    skewness::Float64       # skewness (3rd moment)
    z_statistic::Float64    # z-statistic
end


function agostino_stat(x::AbstractVector{<:Real})
    
    #x = sort(x)
    n = length(x)

    if n < 8 || n > 46340
        throw(ArgumentError("Sample size must be ≥ 8 and ≤ 46340"))
    end
    
    g₁ = skewness(x)
    Y = g₁ * sqrt((n + 1) * (n + 3) / (6 * (n - 2)))
    β₂ = 3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3) / ((n - 2) * (n + 5) * (n + 7) * (n + 9))
    W² = -1 + sqrt(2 * (β₂ - 1))
    δ = 1 / sqrt(log(sqrt(W²)))
    α = sqrt(2 / (W² - 1))
    Z = δ * log(Y / α + sqrt((Y / α)^2 + 1))

    return (n, g₁, Z)    
end


"""
    AgostinoTest(x::AbstractVector{<:Real})

Performs D'Agostino's test of the null hypothesis that the data in vector `x` are normally distributed, against the alternative that the data are skewed. The test is only applicable for N ≥ 8, and large sample sizes (e.g. N > 46340) can result in √-1 in the calculations. Hence, the length of `x` is restricted to be between 8 and 46340 (inclusive). `x` should not contain any missing values.

# References
D'Agostino, R.B. (1970). Transformation to Normality of the Null Distribution of g_1. Biometrika, 57(3) 679-681.

Implements: [`pvalue`](@ref)
"""
function AgostinoTest(x::AbstractVector{<:Real})
    AgostinoTest(agostino_stat(x)...)    
end


testname(::AgostinoTest) = "D'Agostino's test for skewness"

# parameter of interest: name, value under h0, point estimate
population_param_of_interest(x::AgostinoTest) = ("Skewness", 0.0, x.skewness)

function show_params(io::IO, x::AgostinoTest, ident = "")
    println(io, ident, "number of observations: $(x.nobs)")
    println(io, ident, "z-statistic:            $(round(x.z_statistic, digits = 3))")
end

StatsAPI.nobs(x::AgostinoTest) = x.nobs

function StatsAPI.pvalue(x::AgostinoTest; tail = :both)

    check_tail(tail)
    
    left_tail = cdf(Normal(), x.z_statistic)
    right_tail = ccdf(Normal(), x.z_statistic)
    
    if tail == :both
        p = minimum(2 * [left_tail, right_tail])
    elseif tail == :left
        p = left_tail
    else tail == :right
        p = right_tail
    end

    return p
end

