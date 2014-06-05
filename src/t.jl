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

# FIXME(cs): https://github.com/JuliaLang/julia/issues/4935
abstract TTest <: HypothesisTest
abstract TwoSampleTTest <: TTest

pvalue(x::TTest; tail=:both) = pvalue(TDist(x.df), x.t; tail=tail)

# confidence interval by inversion
function ci(x::TTest, alpha::Float64=0.05; tail=:both)
    check_alpha(alpha)

    if tail == :left
        (-Inf, ci(x, alpha*2)[2])
    elseif tail == :right
        (ci(x, alpha*2)[1], Inf)
    elseif tail == :both
        q = quantile(TDist(x.df), 1-alpha/2)
        (x.xbar-q*x.stderr, x.xbar+q*x.stderr)
    else
        error("tail=$(tail) is invalid")
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

function OneSampleTTest(xbar::Real, stddev::Real, n::Int, μ0::Real=0)
    stderr = stddev/sqrt(n)
    t = (xbar-μ0)/stderr
    df = n-1
    OneSampleTTest(n, xbar, df, stderr, t, μ0)
end

OneSampleTTest{T<:Real}(v::AbstractVector{T}, μ0::Real=0) = OneSampleTTest(mean(v), std(v), length(v), μ0)

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

function UnequalVarianceTTest{T<:Real,S<:Real}(x::AbstractVector{T}, y::AbstractVector{S}, μ0::Real=0)
    nx, ny = length(x), length(y)
    xbar = mean(x)-mean(y)
    varx, vary = var(x), var(y)
    stderr = sqrt(varx/nx + vary/ny)
    t = (xbar-μ0)/stderr
    df = (varx / nx + vary / ny)^2 / ((varx / nx)^2 / (nx - 1) + (vary / ny)^2 / (ny - 1))
    UnequalVarianceTTest(nx, ny, xbar, df, stderr, t, μ0)
end

