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

@compat abstract type TTest <: HypothesisTest end
@compat abstract type TwoSampleTTest <: TTest end

pvalue(x::TTest; tail=:both) = pvalue(TDist(x.df), x.t; tail=tail)

default_tail(test::TTest) = :both

# confidence interval by inversion
"""
    confint(test::HypothesisTest, alpha = 0.05; tail = :both)

Compute a confidence interval C with coverage 1-`alpha`.

If `tail` is `:both` (default), then a two-sided confidence interval is returned. If `tail`
is `:left` or `:right`, then a one-sided confidence interval is returned.

!!! note
    Most of the implemented confidence intervals are *strongly consistent*, that is, the
    confidence interval with coverage 1-`alpha` does not contain the test statistic under
    ``h_0`` if and only if the corresponding test rejects the null hypothesis
    ``h_0: \\theta=\\theta_0``:
    ```math
        C (x, 1 − \\alpha) = \\{\\theta : p_\\theta (x) > \\alpha\\},
    ```
    where ``p_\\theta`` is the [`pvalue`](@ref) of the corresponding test.
"""
function StatsBase.confint(x::TTest, alpha::Float64=0.05; tail=:both)
    check_alpha(alpha)

    if tail == :left
        (-Inf, StatsBase.confint(x, alpha*2)[2])
    elseif tail == :right
        (StatsBase.confint(x, alpha*2)[1], Inf)
    elseif tail == :both
        q = quantile(TDist(x.df), 1-alpha/2)
        (x.xbar-q*x.stderr, x.xbar+q*x.stderr)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end


## ONE SAMPLE T-TEST

immutable OneSampleTTest <: TTest
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
    OneSampleTTest(xbar::Real, stdev::Real, n::Int, mu0::Real = 0)

Perform a one sample t-test of the null hypothesis that `n` values with mean `xbar` and
sample standard deviation `stdev`  come from a distribution with `mu0` against the
alternative hypothesis that the distribution does not have mean `mu0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function OneSampleTTest(xbar::Real, stddev::Real, n::Int, μ0::Real=0)
    stderr = stddev/sqrt(n)
    t = (xbar-μ0)/stderr
    df = n-1
    OneSampleTTest(n, xbar, df, stderr, t, μ0)
end

"""
    OneSampleTTest(v::AbstractVector{T<:Real}, mu0::Real = 0)

Perform a one sample t-test of the null hypothesis that the data in vector `v` comes from
a distribution with mean `mu0` against the alternative hypothesis that the distribution
does not have mean `mu0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
OneSampleTTest{T<:Real}(v::AbstractVector{T}, μ0::Real=0) = OneSampleTTest(mean(v), std(v), length(v), μ0)

"""
    OneSampleTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, mu0::Real = 0)

Perform a paired sample t-test of the null hypothesis that the differences between pairs of
values in vectors `x` and `y` come from a distribution with mean `mu0` against the
alternative hypothesis that the distribution does not have mean `mu0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function OneSampleTTest{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0)
    check_same_length(x, y)

    OneSampleTTest(x - y, μ0)
end


## TWO SAMPLE T-TEST (EQUAL VARIANCE)

immutable EqualVarianceTTest <: TwoSampleTTest
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
    EqualVarianceTTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

Perform a two-sample t-test of the null hypothesis that `x` and `y` come from a
distributions with the same mean and equal variances against the alternative hypothesis
that the distributions have different means and but equal variances.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function EqualVarianceTTest{T<:Real,S<:Real}(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0)
    nx, ny = length(x), length(y)
    xbar = mean(x) - mean(y)
    stddev = sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
    stderr = stddev * sqrt(1/nx + 1/ny)
    t = (xbar - μ0) / stderr
    df = nx + ny - 2
    EqualVarianceTTest(nx, ny, xbar, df, stderr, t, μ0)
end


## TWO SAMPLE T-TEST (UNEQUAL VARIANCE)

immutable UnequalVarianceTTest <: TwoSampleTTest
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
from a distributions with the same mean against the alternative hypothesis that the
distributions have different means.

This test is also known as sometimes known as Welch's t-test. It differs from the equal
variance t-test in that it computes the number of degrees of freedom of the test using the
Welch-Satterthwaite equation:
```math
    \\nu_{\\chi'} \\approx \\frac{\\left(\\sum_{i=1}^n k_i s_i^2\\right)^2}{\\sum_{i=1}^n
        \\frac{(k_i s_i^2)^2}{\\nu_i}}
```

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function UnequalVarianceTTest{T<:Real,S<:Real}(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0)
    nx, ny = length(x), length(y)
    xbar = mean(x)-mean(y)
    varx, vary = var(x), var(y)
    stderr = sqrt(varx/nx + vary/ny)
    t = (xbar-μ0)/stderr
    df = (varx / nx + vary / ny)^2 / ((varx / nx)^2 / (nx - 1) + (vary / ny)^2 / (ny - 1))
    UnequalVarianceTTest(nx, ny, xbar, df, stderr, t, μ0)
end
