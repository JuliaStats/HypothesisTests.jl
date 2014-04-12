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

abstract TTest <: HypothesisTest
abstract TwoSampleTTest <: TTest

pvalue(x::TTest; tail=:both) = pvalue(TDist(x.df), x.t; tail=tail)

function ci(x::TTest, alpha::Float64=0.05; tail=:both)
    check_alpha(alpha)
    if tail == :both
        alpha /= 2
    end

    q = quantile(TDist(x.df), alpha)

    if tail == :both
        ((x.t+q)*x.s, (x.t-q)*x.s)
    elseif tail == :left
        (-Inf, (x.t-q)*x.s)
    elseif tail == :right
        ((x.t+q)*x.s, Inf)
    else
        error("tail=$(tail) is invalid")
    end
end

## ONE SAMPLE T-TEST

immutable OneSampleTTest{T<:Real} <: TTest
    t::Float64
    df::T
    s::Float64
end
function OneSampleTTest(xbar::Real, stdev::Real, n::Int, mu0::Real=0)
    s = stdev/sqrt(n)
    OneSampleTTest((xbar-mu0)/s, n-1, s)
end
OneSampleTTest{T<:Real}(v::AbstractVector{T}, mu0::Real=0) =
    OneSampleTTest(mean(v), std(v), length(v)) 
function OneSampleTTest{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S}, mu0::Real=0)
    check_same_length(x, y)
    OneSampleTTest(x - y, mu0)
end
testname(::OneSampleTTest) = "One sample t-test"

## EQUAL VARIANCE T-TEST

immutable EqualVarianceTTest <: TwoSampleTTest
    t::Float64
    df::Int
    s::Float64
end
function EqualVarianceTTest{T<:Real,S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})
    s = sqrt(((length(x)-1) * var(x) + (length(y)-1) * var(y)) /
        (length(x)+length(y)-2) * (1/length(x)+1/length(y)))
    EqualVarianceTTest((mean(x) - mean(y))/s, length(x) + length(y) - 2, s)
end

testname(::EqualVarianceTTest) = "Equal variance t-test"

## UNEQUAL VARIANCE T-TEST

immutable UnequalVarianceTTest <: TwoSampleTTest
    t::Float64
    df::Float64
    s::Float64
end
function UnequalVarianceTTest{T<:Real,S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})
    nx = length(x)
    ny = length(y)
    varx = var(x)
    vary = var(y)
    
    s = sqrt(varx/nx + vary/ny)
    UnequalVarianceTTest((mean(x) - mean(y))/s,
        varx == vary == 0
        ? nx + ny - 2
        : (varx/nx + vary/ny)^2/((varx/nx)^2/(nx - 1) + (vary/ny)^2/(ny - 1)),
        s)
end

testname(::UnequalVarianceTTest) = "Unequal variance t-test"
