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

export OneSampleTTest, EqualVarianceTTest, UnequalVarianceTTest, df

abstract TTest

## ONE SAMPLE T-TEST

type OneSampleTTest{T <: Real} <: TTest
    t::Float64
    df::T
    p_value::Float64
end
test_name(::Type{OneSampleTTest}) = "One sample t-test"

# Based on t and degrees of freedom
function p_value(::Type{OneSampleTTest}, t::Real, df::Real)
    dist = TDist(df)
    2 * min(ccdf(dist, t), cdf(dist, t))
end
left_p_value(::Type{OneSampleTTest}, t::Real, df::Real) = cdf(TDist(df), t)
right_p_value(::Type{OneSampleTTest}, t::Real, df::Real) = ccdf(TDist(df), t)
OneSampleTTest(t::Real, df::Real) = OneSampleTTest(t, df, p_value(OneSampleTTest, t, df))

# Alternative forms
test_statistic{T <: Real}(::Type{OneSampleTTest}, v::Vector{T}, mu::Real) =
    (mean(v)-mu)/(std(v)/sqrt(length(v)))
test_statistic{T <: Real}(::Type{OneSampleTTest}, v::Vector{T}) = test_statistic(OneSampleTTest, v, 0)
function test_statistic{T <: Real, S <: Real}(::Type{OneSampleTTest}, x::Vector{T}, y::Vector{S}, mu::Real)
    if length(x) != length(y)
        error("x and y must be the same length")
    end
    test_statistic(OneSampleTTest, x - y, mu)
end
test_statistic{T <: Real, S <: Real}(::Type{OneSampleTTest}, x::Vector{T}, y::Vector{S}) =
    test_statistic(OneSampleTTest, x, y, 0)

df{T <: Real}(v::Vector{T}) = length(v) - 1
function df{T <: Real, S <: Real}(x::Vector{T}, y::Vector{S})
    if length(x) != length(y)
        error("x and y must be the same length")
    end
    length(x) -1
end

for fn in (:p_value, :left_p_value, :right_p_value)
    @eval begin
        $(fn){T <: Real}(::Type{OneSampleTTest}, v::Vector{T}, varargs...) =
            $(fn)(OneSampleTTest, test_statistic(OneSampleTTest, v, varargs...), length(v) - 1)
    end
end

function OneSampleTTest{T <: Real}(v::Vector{T}, varargs...)
    t = ftest_statistic(OneSampleTTest, v, varargs...)
    df = length(v) - 1
    OneSampleTTest(t, df, p_value(OneSampleTTest, t, df))
end

## EQUAL VARIANCE T-TEST

abstract TwoSampleTTest <: TTest
type EqualVarianceTTest <: TwoSampleTTest
    sd::Float64
    t::Float64
    df::Int
    p_value::Float64
end

sd{T <: Real, S <: Real}(::Type{EqualVarianceTTest}, x::Vector{T}, y::Vector{S}) =
    sqrt(((length(x)-1)*var(x)+(length(y)-1)*var(y))/(length(x)+length(y)-2))

test_statistic{T <: Real, S <: Real}(::Type{EqualVarianceTTest}, x::Vector{T}, y::Vector{S}, sd::Real) =
    (mean(x) - mean(y))/(sd*sqrt(1/length(x)+1/length(y)))
test_statistic{T <: Real, S <: Real}(::Type{EqualVarianceTTest}, x::Vector{T}, y::Vector{S}) =
    test_statistic(EqualVarianceTTest, x, y, sd(EqualVarianceTTest, x, y))

df{T <: Real, S <: Real}(::Type{EqualVarianceTTest}, x::Vector{T}, y::Vector{S}) = length(x) + length(y) - 2

for fn in (:p_value, :left_p_value, :right_p_value)
    @eval begin
        $(fn){T <: Real, S <: Real}(::Type{EqualVarianceTTest}, x::Vector{T}, y::Vector{S}) =
            $(fn)(OneSampleTTest, test_statistic(EqualVarianceTTest, x, y), length(x) + length(y) - 2)
    end
end

function EqualVarianceTTest{T <: Real, S <: Real}(x::Vector{T}, y::Vector{S})
    std = sd(EqualVarianceTTest, x, y)
    t = test_statistic(EqualVarianceTTest, x, y, std)
    dof = length(x) + length(y) - 2
    EqualVarianceTTest(std, t, dof, p_value(OneSampleTTest, t, dof))
end

## UNEQUAL VARIANCE T-TEST

type UnequalVarianceTTest <: TwoSampleTTest
    sd::Float64
    t::Float64
    df::Float64
    p_value::Float64
end

sd(::Type{UnequalVarianceTTest}, nx::Int, varx::Real, ny::Int, vary::Real) =
    sqrt(varx/nx + vary/ny)

test_statistic{T <: Real, S <: Real}(::Type{UnequalVarianceTTest}, x::Vector{T}, y::Vector{S}) =
    (mean(x) - mean(y))/sd(UnequalVarianceTTest, length(x), var(x), length(y), var(y))

df(::Type{UnequalVarianceTTest}, nx::Int, varx::Real, ny::Int, vary::Real) = 
    (varx/nx + vary/ny)^2/((varx/nx)^2/(nx - 1) + (vary/ny)^2/(ny - 1))
df{T <: Real, S <: Real}(::Type{UnequalVarianceTTest}, x::Vector{T}, y::Vector{S}) =
    df(UnequalVarianceTTest, length(x), var(x), length(y), var(y))

for fn in (:p_value, :left_p_value, :right_p_value)
    @eval begin
        function $(fn){T <: Real, S <: Real}(::Type{UnequalVarianceTTest}, x::Vector{T}, y::Vector{S})
            args = (length(x), var(x), length(y), var(y))
            $(fn)(OneSampleTTest, (mean(x) - mean(y))/sd(UnequalVarianceTTest, args...), df(UnequalVarianceTTest, args...))
        end
    end
end

function UnequalVarianceTTest{T <: Real, S <: Real}(x::Vector{T}, y::Vector{S})
    args = (length(x), var(x), length(y), var(y))
    std = sd(UnequalVarianceTTest, args...)
    t = (mean(x) - mean(y))/std
    dof = df(UnequalVarianceTTest, args...)
    UnequalVarianceTTest(std, t, dof, p_value(OneSampleTTest, t, dof))
end