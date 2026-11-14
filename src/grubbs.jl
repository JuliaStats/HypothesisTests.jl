# Grubbs test for detecting outliers
export GrubbsTest

using StatsBase: zscore
using Distributions: TDist, ccdf

struct GrubbsTest{R<:Real} <: HypothesisTest
    n::Int
    G::R
    i::Int
    x_i::R
    tail::Symbol
end

testname(::GrubbsTest) = "Grubbs' test for outliers"

function population_param_of_interest(t::GrubbsTest)
    param, explanation = if t.tail == :both
        "Maximum absolute z-score (G)", "large or small"
    elseif t.tail == :right
        "Maximum z-score (G)", "large"
    else
        "Minimum z-score (G)", "small"
    end
    return param, "All data come from a normal distribution without any unusually $explanation values.", t.G
end

default_tail(::GrubbsTest) = :both

function show_params(io::IO, t::GrubbsTest, ident)
    println(io, ident, "number of observations: ", t.n)
    println(io, ident, "G-statistic:            ", t.G)
    println(io, ident, "potential outlier at:   index $(t.i) (value $(t.x_i))")
end

function StatsAPI.pvalue(gt::GrubbsTest)
    n, G = gt.n, gt.G
    gt.tail == :left && (G = -G)
    s = G^2 * n * (2 - n) / (G^2 * n - (n - 1)^2)
    p = n * ccdf(TDist(n - 2), sqrt(s))
    p > 1 && (p = one(G))
    if gt.tail == :both
        p *= 2
        p > 1 && (p = 2 - p)
    end
    return p
end

function GrubbsTest(x::AbstractVector{<:Real}; tail::Symbol=:both)
    check_tail(tail)
    n = length(x)
    if n < 3
        throw(ArgumentError("Grubbs' test requires at least 3 observations."))
    end

    z = zscore(x)

    G, i = if tail == :both
        findmax(abs.(z))
    elseif tail == :right
        findmax(z)
    else # :left
        findmin(z)
    end

    GrubbsTest(n, G, i, x[i], tail)
end
