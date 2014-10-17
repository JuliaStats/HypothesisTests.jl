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

abstract KSTest <: HypothesisTest
abstract ApproximateKSTest <: KSTest
abstract ExactKSTest <: KSTest

population_param_of_interest(x::KSTest) = ("Supremum of CDF differences", 0.0, x.δ) # parameter of interest: name, value under h0, point estimate

## ONE SAMPLE KS-TEST

# compute supremum of differences between target and empirical cdf before and after the jump of the empirical cdf.
function ksstats{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution)
    n = length(x)
    cdfs = cdf(d, sort(x))
    δp = maximum([1:n] / n - cdfs)
    δn = -minimum([0:n-1] / n - cdfs)
    δ = max(δn, δp)
    (n, δ, δp, δn)
end

### EXACT KOLMOGOROV SMIRNOV TEST

immutable ExactOneSampleKSTest <: ExactKSTest
    n::Int      # number of observations
    δ::Float64  # supremum of CDF differences
    δp::Float64 # supremum of the positive CDF differences
    δn::Float64 # supremum of the negative CDF differences
end

function ExactOneSampleKSTest{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution)
    if length(x) > length(unique(x))
        warn("This test is inaccurate with ties")
    end

    ExactOneSampleKSTest(ksstats(x, d)...)
end

testname(::ExactOneSampleKSTest) = "Exact one sample Kolmorov-Smirnov test"

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
        error("tail=$(tail) is invalid")
    end
end

### APPROXIMATE KOLMOGOROV SMIRNOV TEST

immutable ApproximateOneSampleKSTest <: ApproximateKSTest
    n::Int      # number of observations
    δ::Float64  # supremum of CDF differences
    δp::Float64 # supremum of the positive CDF differences
    δn::Float64 # suproemum of the negative CDF differences
end

function ApproximateOneSampleKSTest{T<:Real}(x::AbstractVector{T}, d::UnivariateDistribution) 
    if length(x) > length(unique(x))
        warn("This test is inaccurate with ties")
    end

    ApproximateOneSampleKSTest(ksstats(x, d)...)
end

testname(::ApproximateOneSampleKSTest) = "Approximate one sample Kolmorov-Smirnov test"

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
        error("tail=$(tail) is invalid")
    end
end

## TWO SAMPLE KS-TEST

### APPROXIMATE KOLMOGOROV SMIRNOV TEST

immutable ApproximateTwoSampleKSTest <: ApproximateKSTest
    n_x::Int    # number of observations
    n_y::Int    # number of observations
    δ::Float64  # supremum of CDF differences
    δp::Float64 # supremum of the positive CDF differences
    δn::Float64 # suproemum of the negative CDF differences
end

function ApproximateTwoSampleKSTest{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})
    n_x, n_y = length(x), length(y)
    if n_x+n_y > length(unique([x,y]))
        warn("This test is inaccurate with ties")
    end

    ApproximateTwoSampleKSTest(ksstats(x, y)...)
end

testname(::ApproximateTwoSampleKSTest) = "Approximate two sample Kolmorov-Smirnov test"

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
        error("tail=$(tail) is invalid")
    end
end

# compute supremum of differences between empirical cdfs.
function ksstats{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})
    n_x, n_y = length(x), length(y)
    sort_idx = sortperm([x, y])
    pdf_diffs = [ones(n_x)/n_x, -ones(n_y)/n_y][sort_idx]
    cdf_diffs = cumsum(pdf_diffs)
    δp = maximum(cdf_diffs)
    δn = -minimum(cdf_diffs)
    δ = max(δp, δn)
    (n_x, n_y, δ, δp, δn)
end
