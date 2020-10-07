# Correlation
# TODO: Implement `CorrelationTest` for two arguments

using Statistics: clampcor

export CorrelationTest

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
struct CorrelationTest{T<:Real} <: HypothesisTest
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
    string("Test for nonzero ", p.k != 0 ? "partial " : "", "correlation")

function population_param_of_interest(p::CorrelationTest)
    param = p.k != 0 ? "Partial correlation" : "Correlation"
    (param, zero(p.r), p.r)
end

StatsBase.nobs(p::CorrelationTest) = p.n
StatsBase.dof(p::CorrelationTest) = p.n - 2 - p.k

function StatsBase.confint(test::CorrelationTest{T}, level::Float64=0.95) where T
    dof(test) > 1 || return (-one(T), one(T))  # Otherwise we can get NaNs
    q = quantile(Normal(), 1 - (1-level) / 2)
    fisher = atanh(test.r)
    bound = q / sqrt(dof(test) - 1)
    elo = clampcor(tanh(fisher - bound))
    ehi = clampcor(tanh(fisher + bound))
    return (elo, ehi)
end

default_tail(::CorrelationTest) = :both
pvalue(test::CorrelationTest; tail=:both) = pvalue(TDist(dof(test)), test.t, tail=tail)

function show_params(io::IO, test::CorrelationTest, indent="")
    println(io, indent, "number of observations:          ", nobs(test))
    println(io, indent, "number of conditional variables: ", test.k)
    println(io, indent, "t-statistic:                     ", test.t)
    println(io, indent, "degrees of freedom:              ", dof(test))
end
