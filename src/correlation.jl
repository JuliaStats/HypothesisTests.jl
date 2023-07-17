# Correlation
# TODO: Implement `CorrelationTest` for two arguments

using Statistics: clampcor
using StatsBase: corspearman

export CorrelationTest, SpearmanCorrelationTest

abstract type AbstractCorrelationTest <: HypothesisTest end

"""
    CorrelationTest(x, y)

Perform a t-test for the hypothesis that ``\\text{Cor}(x,y) = 0``, i.e. the correlation 
of vectors `x` and `y` is zero.

    CorrelationTest(x, y, Z)

Perform a t-test for the hypothesis that ``\\text{Cor}(x,y|Z=z) = 0``, i.e. the partial
correlation of vectors `x` and `y` given the matrix `Z` is zero.

Implements `pvalue` for the t-test.

Implements `confint` using an approximate confidence interval based on Fisher's
``z``-transform.

See also `partialcor` from StatsBase.

# External resources

* [Partial correlation on Wikipedia](https://en.wikipedia.org/wiki/Partial_correlation#As_conditional_independence_test)
  (for the construction of the confidence interval)
* [Section testing using Student's t-distribution from Pearson correlation coefficient on Wikipedia](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Testing_using_Student's_t-distribution)
* [K.J. Levy and S.C. Narula (1978): Testing Hypotheses concerning Partial Correlations: Some Methods and Discussion. International Statistical Review 46(2).](https://www.jstor.org/stable/1402814)
"""
struct CorrelationTest{T<:Real} <: AbstractCorrelationTest
    r::T
    n::Int
    k::Int
    t::T

    # Error checking is done in `cor`
    function CorrelationTest(x::AbstractVector, y::AbstractVector)
        r = cor(x, y)
        n = length(x)
        t = r * sqrt((n - 2) / (1 - r^2))
        return new{typeof(r)}(r, n, 0, t)
    end

    # Error checking is done in `partialcor`
    function CorrelationTest(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix)
        r = partialcor(x, y, Z)
        n, k = size(Z)
        t = r * sqrt((n - 2 - k) / (1 - r^2))
        return new{typeof(r)}(r, n, k, t)
    end

    function CorrelationTest(x::AbstractVector, y::AbstractVector, z::AbstractVector)
        r = partialcor(x, y, z)
        n = length(z)
        t = r * sqrt((n - 3) / (1 - r^2))
        return new{typeof(r)}(r, n, 1, t)
    end
end

testname(p::CorrelationTest) =
    string("Test for nonzero ", p.k != 0 ? "partial " : "", "Pearson correlation")

function population_param_of_interest(p::CorrelationTest)
    param = p.k != 0 ? "Partial Pearson correlation" : "Pearson correlation"
    (param, zero(p.r), p.r)
end

StatsAPI.nobs(p::AbstractCorrelationTest) = p.n
StatsAPI.dof(p::AbstractCorrelationTest) = p.n - 2 - p.k

function StatsAPI.confint(test::CorrelationTest{T}, level::Float64=0.95) where T
    dof(test) > 1 || return (-one(T), one(T))  # Otherwise we can get NaNs
    q = quantile(Normal(), 1 - (1-level) / 2)
    fisher = atanh(test.r)
    bound = q / sqrt(dof(test) - 1)
    elo = clampcor(tanh(fisher - bound))
    ehi = clampcor(tanh(fisher + bound))
    return (elo, ehi)
end

default_tail(::AbstractCorrelationTest) = :both
StatsAPI.pvalue(test::AbstractCorrelationTest; tail=:both) = pvalue(TDist(dof(test)), test.t, tail=tail)

function show_params(io::IO, test::AbstractCorrelationTest, indent="")
    println(io, indent, "number of observations:          ", nobs(test))
    println(io, indent, "number of conditional variables: ", test.k)
    println(io, indent, "t-statistic:                     ", test.t)
    println(io, indent, "degrees of freedom:              ", dof(test))
end

"""
    SpearmanCorrelationTest(x, y)

Perform a t-test for the hypothesis that ``\\text{Cor}(x,y) = 0``, i.e. the rank-based Spearman correlation 
of vectors `x` and `y` is zero.

Implements `pvalue` for the t-test.

Implements `confint` using an approximate confidence interval adjusting for the non-normality of the ranks based on [1]. This is still an approximation and which performs insufficient in the case of:

* small sample sizes n < 25
* a high true population Spearman correlation

In these cases a bootstrap confidence interval can perform better [2].

# External resources
[1] D. G. Bonett and T. A. Wright, “Sample size requirements for estimating pearson, kendall and spearman correlations,” Psychometrika, vol. 65, no. 1, pp. 23–28, Mar. 2000, doi: 10.1007/BF02294183.

[2] A. J. Bishara and J. B. Hittner, “Confidence intervals for correlations when data are not normal,” Behav Res, vol. 49, no. 1, pp. 294–309, Feb. 2017, doi: 10.3758/s13428-016-0702-8.

"""
struct SpearmanCorrelationTest{T<:Real} <: AbstractCorrelationTest
    r::T
    n::Int
    k::Int
    t::T

    # Error checking is done in `corspearman`
    function SpearmanCorrelationTest(x::AbstractVector, y::AbstractVector)
        r = corspearman(x, y)
        n = length(x)
        t = r * sqrt((n - 2) / (1 - r^2))
        return new{typeof(r)}(r, n, 0, t)
    end
end

testname(p::SpearmanCorrelationTest) =  "Spearman correlation"

function population_param_of_interest(p::SpearmanCorrelationTest)
    param = "Spearman Correlation"
    (param, zero(p.r), p.r)
end

function StatsAPI.confint(test::SpearmanCorrelationTest{T}, level::Float64=0.95) where T
    dof(test) > 1 || return (-one(T), one(T))  # Otherwise we can get NaNs
    q = quantile(Normal(), 1 - (1-level) / 2)
    fisher = atanh(test.r)
    bound = sqrt((1 + test.r^2 / 2) / (dof(test)-1)) * q # Estimates variance as in Bonnet et al. (2000)
    elo = clampcor(tanh(fisher - bound))
    ehi = clampcor(tanh(fisher + bound))
    return (elo, ehi)
end