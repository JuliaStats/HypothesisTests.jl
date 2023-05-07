# binomial.jl
# Binomial tests
#
# Copyright (C) 2013   Simon Kornblith
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

export BinomialTest, SignTest

## BINOMIAL TEST

struct BinomialTest <: HypothesisTest
    p::Float64
    x::Int
    n::Int

    BinomialTest(x::Real, n::Real, p::Real=0.5) = new(p, x, n)
end

"""
    BinomialTest(x::Integer, n::Integer, p::Real = 0.5)
    BinomialTest(x::AbstractVector{Bool}, p::Real = 0.5)

Perform a binomial test of the null hypothesis that the distribution from which `x`
successes were encountered in `n` draws (or alternatively from which the vector `x` was
drawn) has success probability `p` against the alternative hypothesis that the success
probability is not equal to `p`.

Computed confidence intervals by default are Clopper-Pearson intervals.
See the [`confint(::BinomialTest)`](@ref) documentation for a list of
supported methods to compute confidence intervals.

Implements: [`pvalue`](@ref), [`confint(::BinomialTest)`](@ref)
"""
BinomialTest(x::AbstractVector{Bool}, p=0.5) =
    BinomialTest(sum(x), length(x), p)

"""
    testname(::HypothesisTest)

Returns the string value, e.g. "Binomial test" or "Sign Test".
"""
testname(::BinomialTest) = "Binomial test"
population_param_of_interest(x::BinomialTest) = ("Probability of success", x.p, x.x/x.n) # parameter of interest: name, value under h0, point estimate
default_tail(test::BinomialTest) = :both

function show_params(io::IO, x::BinomialTest, ident="")
    println(io, ident, "number of observations: $(x.n)")
    println(io, ident, "number of successes:    $(x.x)")
end

pvalue(x::BinomialTest; tail=:both) = pvalue(Binomial(x.n, x.p), x.x; tail=tail)

# Confidence interval
"""
    confint(test::BinomialTest; level = 0.95, tail = :both, method = :clopper_pearson)

Compute a confidence interval with coverage `level` for a binomial proportion using one
of the following methods. Possible values for `method` are:

  - `:clopper_pearson` (default): Clopper-Pearson interval is based on the binomial
    distribution. The empirical coverage is never less than the nominal coverage of
    `level`; it is usually too conservative.
  - `:wald`: Wald (or normal approximation) interval relies on the standard approximation of
    the actual binomial distribution by a normal distribution. Coverage can be erratically
    poor for success probabilities close to zero or one.
  - `:wilson`: Wilson score interval relies on a normal approximation. In contrast to `:wald`,
    the standard deviation is not approximated by an empirical estimate, resulting in good
    empirical coverages even for small numbers of draws and extreme success probabilities.
  - `:jeffrey`: Jeffreys interval is a Bayesian credible interval obtained by using a
    non-informative Jeffreys prior. The interval is very similar to the Wilson interval.
  - `:agresti_coull`: Agresti-Coull interval is a simplified version of the Wilson interval;
    both are centered around the same value. The Agresti Coull interval has higher or equal
    coverage.
  - `:arcsine`: Confidence interval computed using the arcsine transformation to make
    ``var(p)`` independent of the probability ``p``.

# References

  * Brown, L.D., Cai, T.T., and DasGupta, A. Interval estimation for a binomial proportion.
    Statistical Science, 16(2):101–117, 2001.

# External links

  * [Binomial confidence interval on Wikipedia](https://en.wikipedia.org/wiki/
    Binomial_proportion_confidence_interval)
"""
function StatsBase.confint(x::BinomialTest; level::Float64=0.95, tail=:both, method=:clopper_pearson)
    check_level(level)

    if tail == :left
        (0.0, StatsBase.confint(x, level=1-(1-level)*2, method=method)[2])
    elseif tail == :right
        (StatsBase.confint(x, level=1-(1-level)*2, method=method)[1], 1.0)
    elseif tail == :both
        if method == :clopper_pearson
            ci_clopper_pearson(x, 1-level)
        elseif method == :wald
            ci_wald(x, 1-level)
        elseif method == :wilson
            ci_wilson(x, 1-level)
        elseif method == :jeffrey
            ci_jeffrey(x, 1-level)
        elseif method == :agresti_coull
            ci_agresti_coull(x, 1-level)
        elseif method == :arcsine
            ci_arcsine(x, 1-level)
        else
            throw(ArgumentError("method=$(method) is not implemented yet"))
        end
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

# Clopper-Pearson interval (confidence interval by inversion)
function ci_clopper_pearson(x::BinomialTest, alpha::Float64=0.05)
    (x.x == 0 ? 0.0 : quantile(Beta(x.x, x.n - x.x + 1), alpha/2),
     x.x == x.n ? 1.0 : quantile(Beta(x.x + 1, x.n - x.x), 1-alpha/2))
end

# Wald interval / normal approximation interval
function ci_wald(x::BinomialTest, alpha::Float64=0.05)
    μ = x.x / x.n
    σ = sqrt(μ*(1-μ)/x.n)
    lower, upper = (quantile(Normal(μ, σ), alpha/2), quantile(Normal(μ, σ), 1-alpha/2))
    # make sure we stay in [0, 1]
    (max(lower, 0), min(upper, 1))
end

# Jeffreys interval
function ci_jeffrey(x::BinomialTest, alpha::Float64=0.05)
    (quantile(Beta(x.x + 1/2, x.n - x.x + 1/2), alpha/2), quantile(Beta(x.x + 1/2, x.n - x.x + 1/2), 1-alpha/2))
end

# Agresti-Coull interval
function ci_agresti_coull(x::BinomialTest, alpha::Float64=0.05)
    q = quantile(Normal(), 1-alpha/2)
    n = x.n + q^2
    μ = (x.x + q^2/2)/n
    σ = sqrt(μ*(1-μ)/n)
    # make sure we stay in [0, 1]
    (max(μ-q*σ, 0), min(μ+q*σ, 1))
end

# Wilson score interval
function ci_wilson(x::BinomialTest, alpha::Float64=0.05)
    q = quantile(Normal(), 1-alpha/2)
    p = x.x / x.n
    denominator = 1 + q^2/x.n
    μ = p + q^2/(2*x.n)
    μ /= denominator
    σ = sqrt(p*(1-p)/x.n + q^2/(4x.n^2))
    σ /= denominator
    # make sure we stay in [0, 1]
    (max(μ-q*σ, 0), min(μ+q*σ, 1))
end

# Arcsine transformation interval as based on Cohen's H: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Arcsine_transformation
function ci_arcsine(x::BinomialTest, alpha::Float64=0.05)
    q = quantile(Normal(), 1-alpha/2)
    p = x.x / x.n
    z = q/(2*sqrt(x.n))
    (sin(asin(sqrt(p))-z)^2, sin(asin(sqrt(p))+z)^2)
end

## SIGN TEST

struct SignTest <: HypothesisTest
    median::Float64
    x::Int
    n::Int
    data
end

"""
    SignTest(x::AbstractVector{T<:Real}, median::Real = 0)
    SignTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, median::Real = 0)

Perform a sign test of the null hypothesis that the distribution from which `x`
(or `x - y` if `y` is provided) was drawn has median `median` against the alternative
hypothesis that the median is not equal to `median`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
SignTest(x::AbstractVector{T}, median::Real=0) where {T<:Real} =
    SignTest(median, sum(x .> median), sum(x .!= median), sort(x))
SignTest(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real} =
    SignTest(x - y, 0.0)

testname(::SignTest) = "Sign Test"
population_param_of_interest(x::SignTest) = ("Median", x.median, median(x.data)) # parameter of interest: name, value under h0, point estimate
default_tail(test::SignTest) = :both

function show_params(io::IO, x::SignTest, ident="")
    text1 = "number of observations:"
    text2 = "observations larger than $(x.median): "
    maxlen = length(text2)

    println(io, ident, text1, repeat(" ", maxlen-length(text1)), x.n)
    println(io, ident, text2, x.x)
end

pvalue(x::SignTest; tail=:both) = pvalue(Binomial(x.n, 0.5), x.x; tail=tail)

function StatsBase.confint(x::SignTest; level::Float64=0.95, tail=:both)
    check_level(level)

    if tail == :left
        q = Int(quantile(Binomial(x.n, 0.5), 1-level))
        (x.data[q+1], x.median)
    elseif tail == :right
        q = Int(quantile(Binomial(x.n, 0.5), 1-level))
        (x.median, x.data[end-q])
    elseif tail == :both
        q = Int(quantile(Binomial(x.n, 0.5), (1-level)/2))
        (x.data[q+1], x.data[end-q])
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end
