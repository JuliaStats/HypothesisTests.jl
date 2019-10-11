# Partial correlation

using Statistics: clampcor

export PartialCorTest

"""
    PartialCorTest(x, y, Z)

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
struct PartialCorTest{T<:Real} <: HypothesisTest
    r::T
    n::Int
    k::Int
    t::T

    # Error checking is done in `partialcor`
    function PartialCorTest(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix)
        r = partialcor(x, y, Z)
        n, k = size(Z)
        t = r * sqrt((n - 2 - k) / (1 - r^2))
        return new{typeof(r)}(r, n, k, t)
    end

    function PartialCorTest(x::AbstractVector, y::AbstractVector, z::AbstractVector)
        r = partialcor(x, y, z)
        n = length(z)
        t = r * sqrt((n - 3) / (1 - r^2))
        return new{typeof(r)}(r, n, 1, t)
    end
end

testname(::PartialCorTest) = "Test for partial correlation"
population_param_of_interest(p::PartialCorTest) = ("Partial correlation", zero(p.r), p.r)

StatsBase.nobs(p::PartialCorTest) = p.n
StatsBase.dof(p::PartialCorTest) = p.n - 2 - p.k

function StatsBase.confint(test::PartialCorTest{T}, alpha::Float64=0.05) where T
    dof(test) > 1 || return (-one(T), one(T))  # Otherwise we can get NaNs
    q = quantile(Normal(), 1 - alpha / 2)
    fisher = atanh(test.r)
    bound = q / sqrt(test.n - 3 - test.k)
    elo = clampcor(tanh(fisher - bound))
    ehi = clampcor(tanh(fisher + bound))
    return (elo, ehi)
end

default_tail(::PartialCorTest) = :both
pvalue(test::PartialCorTest; tail=:both) = pvalue(TDist(dof(test)), test.t, tail=tail)

function show_params(io::IO, test::PartialCorTest, indent="")
    println(io, indent, "number of observations:          ", nobs(test))
    println(io, indent, "number of conditional variables: ", test.k)
    println(io, indent, "t-statistic:                     ", test.t)
    println(io, indent, "degrees of freedom:              ", dof(test))
end
