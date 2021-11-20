# t.jl
# Various forms of t-tests
#
# Copyright (C) 2012   Simon Kornblith
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

export OneSampleTTest, TwoSampleTTest, EqualVarianceTTest,
    UnequalVarianceTTest

abstract type TTest <: HypothesisTest end
abstract type TwoSampleTTest <: TTest end

pvalue(x::TTest; tail=:both) = pvalue(TDist(x.df), x.t; tail=tail)

default_tail(test::TTest) = :both

# confidence interval by inversion
function StatsBase.confint(x::TTest; level::Float64=0.95, tail=:both)
    check_level(level)

    if tail == :left
        (-Inf, StatsBase.confint(x, level=1-(1-level)*2)[2])
    elseif tail == :right
        (StatsBase.confint(x, level=1-(1-level)*2)[1], Inf)
    elseif tail == :both
        q = quantile(TDist(x.df), 1-(1-level)/2)
        (x.xbar-q*x.stderr, x.xbar+q*x.stderr)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end


## ONE SAMPLE T-TEST

struct OneSampleTTest <: TTest
    n::Int       # number of observations
    xbar::Real   # estimated mean
    df::Int      # degrees of freedom
    stderr::Real # empirical standard error
    t::Real      # t-statistic
    μ0::Real     # mean under h_0
end

testname(::OneSampleTTest) = "One sample t-test"
population_param_of_interest(x::OneSampleTTest) = ("Mean", x.μ0, x.xbar) # parameter of interest: name, value under h0, point estimate

function show_params(io::IO, x::OneSampleTTest, ident="")
    println(io, ident, "number of observations:   $(x.n)")
    println(io, ident, "t-statistic:              $(x.t)")
    println(io, ident, "degrees of freedom:       $(x.df)")
    println(io, ident, "empirical standard error: $(x.stderr)")
end

"""
    OneSampleTTest(xbar::Real, stddev::Real, n::Int, μ0::Real = 0)

Perform a one sample t-test of the null hypothesis that `n` values with mean `xbar` and
sample standard deviation `stddev`  come from a distribution with mean `μ0` against the
alternative hypothesis that the distribution does not have mean `μ0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function OneSampleTTest(xbar::Real, stddev::Real, n::Int, μ0::Real=0)
    stderr = stddev/sqrt(n)
    t = (xbar-μ0)/stderr
    df = n-1
    OneSampleTTest(n, xbar, df, stderr, t, μ0)
end

"""
    OneSampleTTest(v::AbstractVector{T<:Real}, μ0::Real = 0)

Perform a one sample t-test of the null hypothesis that the data in vector `v` comes from
a distribution with mean `μ0` against the alternative hypothesis that the distribution
does not have mean `μ0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
OneSampleTTest(v::AbstractVector{T}, μ0::Real=0) where {T<:Real} = OneSampleTTest(mean(v), std(v), length(v), μ0)

"""
    OneSampleTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, μ0::Real = 0)

Perform a paired sample t-test of the null hypothesis that the differences between pairs of
values in vectors `x` and `y` come from a distribution with mean `μ0` against the
alternative hypothesis that the distribution does not have mean `μ0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
    
!!! note
    This test is also known as a t-test for paired or dependent samples, see
    [paired difference test](https://en.wikipedia.org/wiki/Paired_difference_test) on Wikipedia.
"""
function OneSampleTTest(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0) where {T<:Real, S<:Real}
    check_same_length(x, y)

    OneSampleTTest(x - y, μ0)
end


## TWO SAMPLE T-TEST (EQUAL VARIANCE)

struct EqualVarianceTTest <: TwoSampleTTest
    n_x::Int     # number of observations
    n_y::Int     # number of observations
    xbar::Real   # estimated mean difference
    df::Int      # degrees of freedom
    stderr::Real # empirical standard error
    t::Real      # t-statistic
    μ0::Real     # mean difference under h_0
end

function show_params(io::IO, x::TwoSampleTTest, ident="")
    println(io, ident, "number of observations:   [$(x.n_x),$(x.n_y)]")
    println(io, ident, "t-statistic:              $(x.t)")
    println(io, ident, "degrees of freedom:       $(x.df)")
    println(io, ident, "empirical standard error: $(x.stderr)")
end

testname(::EqualVarianceTTest) = "Two sample t-test (equal variance)"
population_param_of_interest(x::TwoSampleTTest) = ("Mean difference", x.μ0, x.xbar) # parameter of interest: name, value under h0, point estimate

"""
    EqualVarianceTTest(nx::Int, ny::Int, mx::Real, my::Real, vx::Real, vy::Real, μ0::Real=0)

Perform a two-sample t-test of the null hypothesis that samples `x` and `y` described by the number
of elements `nx` and `ny`, the mean `mx` and `my`, and variance `vx` and `vy` come from distributions
with equals means and variances. The alternative hypothesis is that the distributions have different
means but equal variances.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function EqualVarianceTTest(nx::Int, ny::Int, mx::Real, my::Real, vx::Real, vy::Real, μ0::Real=0)
    xbar = mx - my
    stddev = sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    stderr = stddev * sqrt(1/nx + 1/ny)
    t = (xbar - μ0) / stderr
    df = nx + ny - 2
    EqualVarianceTTest(nx, ny, xbar, df, stderr, t, μ0)
end

"""
    EqualVarianceTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

Perform a two-sample t-test of the null hypothesis that `x` and `y` come from distributions
with equal means and variances against the alternative hypothesis that the distributions
have different means but equal variances.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function EqualVarianceTTest(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0) where {T<:Real,S<:Real}
    nx, ny = length(x), length(y)
    mx, my = mean(x), mean(y)
    vx, vy = var(x), var(y)
    EqualVarianceTTest(nx, ny, mx, my, vx, vy, μ0)
end


## TWO SAMPLE T-TEST (UNEQUAL VARIANCE)

struct UnequalVarianceTTest <: TwoSampleTTest
    n_x::Int     # number of observations
    n_y::Int     # number of observations
    xbar::Real   # estimated mean
    df::Real     # degrees of freedom
    stderr::Real # empirical standard error
    t::Real      # t-statistic
    μ0::Real     # mean under h_0
end

testname(::UnequalVarianceTTest) = "Two sample t-test (unequal variance)"

"""
    UnequalVarianceTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

Perform an unequal variance two-sample t-test of the null hypothesis that `x` and `y` come
from distributions with equal means against the alternative hypothesis that the
distributions have different means.

This test is sometimes known as Welch's t-test. It differs from the equal variance t-test in
that it computes the number of degrees of freedom of the test using the Welch-Satterthwaite
equation:
```math
    ν_{χ'} ≈ \\frac{\\left(\\sum_{i=1}^n k_i s_i^2\\right)^2}{\\sum_{i=1}^n
        \\frac{(k_i s_i^2)^2}{ν_i}}
```

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function UnequalVarianceTTest(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0) where {T<:Real,S<:Real}
    nx, ny = length(x), length(y)
    xbar = mean(x)-mean(y)
    varx, vary = var(x), var(y)
    stderr = sqrt(varx/nx + vary/ny)
    t = (xbar-μ0)/stderr
    df = (varx / nx + vary / ny)^2 / ((varx / nx)^2 / (nx - 1) + (vary / ny)^2 / (ny - 1))
    UnequalVarianceTTest(nx, ny, xbar, df, stderr, t, μ0)
end
