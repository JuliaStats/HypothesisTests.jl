# kolmogorov_smirnov.jl
# Kolmogorov–Smirnov
#
# Copyright (C) 2014   Christoph Sawade
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

export
    ExactOneSampleKSTest,
    ApproximateOneSampleKSTest, ApproximateTwoSampleKSTest

abstract type KSTest <: HypothesisTest end
abstract type ApproximateKSTest <: KSTest end
abstract type ExactKSTest <: KSTest end

population_param_of_interest(x::KSTest) = ("Supremum of CDF differences", 0.0, x.δ) # parameter of interest: name, value under h0, point estimate
default_tail(test::KSTest) = :both

## ONE SAMPLE KS-TEST

# compute supremum of differences between target and empirical cdf before and after the jump of the empirical cdf.
function ksstats(x::AbstractVector{T}, d::UnivariateDistribution) where T<:Real
    n = length(x)
    cdfs = cdf.(Ref(d), sort(x))
    δp = maximum((1:n) / n - cdfs)
    δn = -minimum((0:n-1) / n - cdfs)
    δ = max(δn, δp)
    (n, δ, δp, δn)
end

### EXACT KOLMOGOROV SMIRNOV TEST

struct ExactOneSampleKSTest <: ExactKSTest
    n::Int      # number of observations
    δ::Float64  # supremum of CDF differences
    δp::Float64 # supremum of the positive CDF differences
    δn::Float64 # supremum of the negative CDF differences
end

"""
    ExactOneSampleKSTest(x::AbstractVector{<:Real}, d::UnivariateDistribution)

Perform a one-sample exact Kolmogorov–Smirnov test of the null hypothesis that the data in
vector `x` comes from the distribution `d` against the alternative hypothesis that the
sample is not drawn from `d`.

Implements: [`pvalue`](@ref)
"""
function ExactOneSampleKSTest(x::AbstractVector{T}, d::UnivariateDistribution) where T<:Real
    if length(x) > length(unique(x))
        @warn("This test is inaccurate with ties")
    end

    ExactOneSampleKSTest(ksstats(x, d)...)
end

testname(::ExactOneSampleKSTest) = "Exact one sample Kolmogorov-Smirnov test"

function show_params(io::IO, x::ExactOneSampleKSTest, ident="")
    println(io, ident, "number of observations:   $(x.n)")
end

function pvalue(x::ExactKSTest; tail=:both)
    if tail == :left
        pvalue(KSOneSided(x.n), x.δn; tail=:right)
    elseif tail == :right
        pvalue(KSOneSided(x.n), x.δp; tail=:right)
    elseif tail == :both
        pvalue(KSDist(x.n), x.δ; tail=:right)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

### APPROXIMATE KOLMOGOROV SMIRNOV TEST

struct ApproximateOneSampleKSTest <: ApproximateKSTest
    n::Int      # number of observations
    δ::Float64  # supremum of CDF differences
    δp::Float64 # supremum of the positive CDF differences
    δn::Float64 # suproemum of the negative CDF differences
end

"""
    ApproximateOneSampleKSTest(x::AbstractVector{<:Real}, d::UnivariateDistribution)

Perform an asymptotic one-sample Kolmogorov–Smirnov test of the null hypothesis that the
data in vector `x` comes from the distribution `d` against the alternative hypothesis
that the sample is not drawn from `d`.

Implements: [`pvalue`](@ref)
"""
function ApproximateOneSampleKSTest(x::AbstractVector{T}, d::UnivariateDistribution) where T<:Real
    if length(x) > length(unique(x))
        @warn("This test is inaccurate with ties")
    end

    ApproximateOneSampleKSTest(ksstats(x, d)...)
end

testname(::ApproximateOneSampleKSTest) = "Approximate one sample Kolmogorov-Smirnov test"

function show_params(io::IO, x::ApproximateOneSampleKSTest, ident="")
    println(io, ident, "number of observations:   $(x.n)")
    println(io, ident, "KS-statistic:             $(sqrt(x.n)*x.δ)")
end

# one-sided: http://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test
function pvalue(x::ApproximateOneSampleKSTest; tail=:both)
    if tail == :left
        exp(-2*x.n*x.δn^2)
    elseif tail == :right
        exp(-2*x.n*x.δp^2)
    elseif tail == :both
        pvalue(Kolmogorov(), sqrt(x.n)*x.δ; tail=:right)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

## TWO SAMPLE KS-TEST

### APPROXIMATE KOLMOGOROV SMIRNOV TEST

struct ApproximateTwoSampleKSTest <: ApproximateKSTest
    n_x::Int    # number of observations
    n_y::Int    # number of observations
    δ::Float64  # supremum of CDF differences
    δp::Float64 # supremum of the positive CDF differences
    δn::Float64 # suproemum of the negative CDF differences
end

"""
    ApproximateTwoSampleKSTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Perform an asymptotic two-sample Kolmogorov–Smirnov-test of the null hypothesis that `x`
and `y` are drawn from the same distribution against the alternative hypothesis that they
come from different distributions.

Implements: [`pvalue`](@ref)

# External links

  * [Approximation of one-sided test (Encyclopedia of Mathematics)
    ](https://www.encyclopediaofmath.org/index.php/Kolmogorov-Smirnov_test)
"""
function ApproximateTwoSampleKSTest(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real}
    n_x, n_y = length(x), length(y)

    ApproximateTwoSampleKSTest(ksstats(x, y)...)
end

testname(::ApproximateTwoSampleKSTest) = "Approximate two sample Kolmogorov-Smirnov test"

function show_params(io::IO, x::ApproximateTwoSampleKSTest, ident="")
    n = x.n_x*x.n_y/(x.n_x+x.n_y)
    println(io, ident, "number of observations:   [$(x.n_x),$(x.n_y)]")
    println(io, ident, "KS-statistic:              $(sqrt(n)*x.δ)")
end

function pvalue(x::ApproximateTwoSampleKSTest; tail=:both)
    n = x.n_x*x.n_y/(x.n_x+x.n_y)
    if tail == :left
        exp(-2*n*x.δn^2)
    elseif tail == :right
        exp(-2*n*x.δp^2)
    elseif tail == :both
        pvalue(Kolmogorov(), sqrt(n)*x.δ; tail=:right)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

# compute supremum of differences between empirical cdfs.
function ksstats(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real}
    n_x, n_y = length(x), length(y)
    all_values = [x; y]
    sort_idx = sortperm(all_values)
    δ_y = 1 / n_y
    δ_x = 1 / n_x
    δ = δp = δn = zero(δ_y)

    for i in 1:(n_x + n_y)
        if sort_idx[i] > n_x
            δ -= δ_y
        else
            δ += δ_x
        end

        # only update δp/δn if the value is about to change or at the last step.
        if i == n_x + n_y || all_values[sort_idx[i]] != all_values[sort_idx[i + 1]]
            if δ > δp
                δp = δ
            elseif δ < δn
                δn = δ
            end
        end
    end

    (n_x, n_y, max(δp, -δn), δp, -δn)
end   
