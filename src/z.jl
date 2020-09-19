# z.jl
# Various forms of z-tests
#
# Copyright (C) 2016   John Myles White
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

export OneSampleZTest, TwoSampleZTest, EqualVarianceZTest,
    UnequalVarianceZTest

abstract type ZTest <: HypothesisTest end
abstract type TwoSampleZTest <: ZTest end

pvalue(x::ZTest; tail=:both) = pvalue(Normal(0.0, 1.0), x.z; tail=tail)

default_tail(test::ZTest) = :both

# confidence interval by inversion
function StatsBase.confint(x::ZTest; level::Float64=0.95, tail=:both)
    check_level(level)

    if tail == :left
        (-Inf, StatsBase.confint(x, level=1-(1-level)*2)[2])
    elseif tail == :right
        (StatsBase.confint(x, level=1-(1-level)*2)[1], Inf)
    elseif tail == :both
        q = cquantile(Normal(0.0, 1.0), (1-level)/2)
        (x.xbar-q*x.stderr, x.xbar+q*x.stderr)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end


## ONE SAMPLE Z-TEST

struct OneSampleZTest <: ZTest
    n::Int       # number of observations
    xbar::Real   # estimated mean
    stderr::Real # population standard error
    z::Real      # z-statistic
    μ0::Real     # mean under h_0
end

testname(::OneSampleZTest) = "One sample z-test"
population_param_of_interest(x::OneSampleZTest) = ("Mean", x.μ0, x.xbar) # parameter of interest: name, value under h0, point estimate

function show_params(io::IO, x::OneSampleZTest, ident="")
    println(io, ident, "number of observations:   $(x.n)")
    println(io, ident, "z-statistic:              $(x.z)")
    println(io, ident, "population standard error: $(x.stderr)")
end

"""
    OneSampleZTest(xbar::Real, stddev::Real, n::Int, μ0::Real = 0)

Perform a one sample z-test of the null hypothesis that `n` values with mean `xbar` and
population standard deviation `stddev`  come from a distribution with mean `μ0` against the
alternative hypothesis that the distribution does not have mean `μ0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function OneSampleZTest(xbar::Real, stddev::Real, n::Int, μ0::Real=0)
    stderr = stddev/sqrt(n)
    z = (xbar-μ0)/stderr
    OneSampleZTest(n, xbar, stderr, z, μ0)
end

"""
    OneSampleZTest(v::AbstractVector{T<:Real}, μ0::Real = 0)

Perform a one sample z-test of the null hypothesis that the data in vector `v` comes from
a distribution with mean `μ0` against the alternative hypothesis that the distribution
does not have mean `μ0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
OneSampleZTest(v::AbstractVector{T}, μ0::Real=0) where {T<:Real} = OneSampleZTest(mean(v), std(v), length(v), μ0)

"""
    OneSampleZTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real}, μ0::Real = 0)

Perform a paired sample z-test of the null hypothesis that the differences between pairs of
values in vectors `x` and `y` come from a distribution with mean `μ0` against the
alternative hypothesis that the distribution does not have mean `μ0`.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function OneSampleZTest(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0) where {T<:Real, S<:Real}
    check_same_length(x, y)

    OneSampleZTest(x - y, μ0)
end


## TWO SAMPLE Z-TEST (EQUAL VARIANCE)

struct EqualVarianceZTest <: TwoSampleZTest
    n_x::Int     # number of observations
    n_y::Int     # number of observations
    xbar::Real   # estimated mean difference
    stderr::Real # population standard error
    z::Real      # z-statistic
    μ0::Real     # mean difference under h_0
end

function show_params(io::IO, x::TwoSampleZTest, ident="")
    println(io, ident, "number of observations:   [$(x.n_x),$(x.n_y)]")
    println(io, ident, "z-statistic:              $(x.z)")
    println(io, ident, "population standard error: $(x.stderr)")
end

testname(::EqualVarianceZTest) = "Two sample z-test (equal variance)"
population_param_of_interest(x::TwoSampleZTest) = ("Mean difference", x.μ0, x.xbar) # parameter of interest: name, value under h0, point estimate

"""
    EqualVarianceZTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

Perform a two-sample z-test of the null hypothesis that `x` and `y` come from distributions
with equal means and variances against the alternative hypothesis that the distributions
have different means but equal variances.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function EqualVarianceZTest(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0) where {T<:Real,S<:Real}
    nx, ny = length(x), length(y)
    xbar = mean(x) - mean(y)
    stddev = sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
    stderr = stddev * sqrt(1/nx + 1/ny)
    z = (xbar - μ0) / stderr
    EqualVarianceZTest(nx, ny, xbar, stderr, z, μ0)
end


## TWO SAMPLE Z-TEST (UNEQUAL VARIANCE)

struct UnequalVarianceZTest <: TwoSampleZTest
    n_x::Int     # number of observations
    n_y::Int     # number of observations
    xbar::Real   # estimated mean
    stderr::Real # empirical standard error
    z::Real      # z-statistic
    μ0::Real     # mean under h_0
end

testname(::UnequalVarianceZTest) = "Two sample z-test (unequal variance)"

"""
    UnequalVarianceZTest(x::AbstractVector{T<:Real}, y::AbstractVector{T<:Real})

Perform an unequal variance two-sample z-test of the null hypothesis that `x` and `y` come
from distributions with equal means against the alternative hypothesis that the
distributions have different means.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function UnequalVarianceZTest(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0) where {T<:Real,S<:Real}
    nx, ny = length(x), length(y)
    xbar = mean(x)-mean(y)
    varx, vary = var(x), var(y)
    stderr = sqrt(varx/nx + vary/ny)
    z = (xbar-μ0)/stderr
    UnequalVarianceZTest(nx, ny, xbar, stderr, z, μ0)
end
