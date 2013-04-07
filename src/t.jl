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

function ci(x::TTest, alpha::Float64; tail=:both)
    if alpha <= 0 || alpha >= 1
        error("alpha $alpha not in range (0, 1)")
    end

    alpha = min(alpha, 1-alpha)

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

immutable OneSampleTTest{T <: Real, S <: Real} <: TTest
    t::Float64
    df::T
    mu0::S
    s::Float64
end
function OneSampleTTest{T <: Real}(v::Vector{T}, mu0::Real)
    s = std(v)/sqrt(length(v))
    OneSampleTTest((mean(v)-mu0)/s, length(v)-1, mu0, s)
end
function OneSampleTTest{T <: Real, S <: Real}(x::Vector{T}, y::Vector{S},
                                              mu0::Real)
    if length(x) != length(y)
        error("x and y must be the same length")
    end
    OneSampleTTest(x - y, mu0)
end
OneSampleTTest{T <: Real}(t::Float64, df::T) = OneSampleTTest(t, df, 0, 1)
OneSampleTTest{T <: Real}(v::Vector{T}) = OneSampleTTest(v, 0)
OneSampleTTest{T <: Real, S <: Real}(x::Vector{T}, y::Vector{S}) =
    OneSampleTTest(x, y, 0)
testname(::OneSampleTTest) = "One sample t-test"

## EQUAL VARIANCE T-TEST

immutable EqualVarianceTTest <: TwoSampleTTest
    t::Float64
    df::Int
    s::Float64
end
function EqualVarianceTTest{T <: Real, S <: Real}(x::Vector{T}, y::Vector{S})
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
function UnequalVarianceTTest{T <: Real, S <: Real}(x::Vector{T}, y::Vector{S})
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

testname(::EqualVarianceTTest) = "Unequal variance t-test"