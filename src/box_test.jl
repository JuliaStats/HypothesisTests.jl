# box_test.jl
# Box-Pierce and Ljung-Box tests for autocorrelation
#
# Copyright (C) 2017   Benjamin Born
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

export BoxPierceTest, LjungBoxTest

# Box-Pierce test

struct BoxPierceTest <: HypothesisTest
    y::Vector{Float64}  # time series vector
    n::Int              # number of observations
    lag::Int            # number of lags in test statistic
    dof::Int            # degrees-of-freedom correction (if y is residual)
    Q::Float64          # test statistic
end

"""
    BoxPierceTest(y, lag, dof=0)

Compute the Box-Pierce `Q` statistic to test the null hypothesis of
independence in a time series `y`.

`lag` specifies the number of lags used in the construction of `Q`. When
testing the residuals of an estimated model, `dof` has to be set to the number
of estimated parameters. E.g., when testing the residuals of an ARIMA(p,0,q)
model, set `dof=p+q`.

# External links

  * [Box-Pierce test on Wikipedia](https://en.wikipedia.org/wiki/Ljung–Box_test#Box-Pierce_test)
"""
function BoxPierceTest(y::AbstractVector{T}, lag::Int, dof::Int=0) where T<:Real
    if dof>=lag
        throw(ArgumentError("Number of lags has to be larger than degrees of" *
        " freedom correction"))
    end
    n = size(y,1)
    Q = n * first(sum(k -> autocor(y,k:k).^2, 1:lag))
    BoxPierceTest(y,n,lag,dof,Q)
end

testname(::BoxPierceTest) = "Box-Pierce autocorrelation test"
population_param_of_interest(x::BoxPierceTest) = ("autocorrelations up to lag k",
    "all zero", NaN)
default_tail(test::BoxPierceTest) = :right

function show_params(io::IO, x::BoxPierceTest, ident)
    println(io, ident, "number of observations:         ", x.n)
    println(io, ident, "number of lags:                 ", x.lag)
    println(io, ident, "degrees of freedom correction:  ", x.dof)
    println(io, ident, "Q statistic:                    ", x.Q)
end

pvalue(x::BoxPierceTest) = pvalue(Chisq(x.lag-x.dof), x.Q; tail=:right)

#Ljung-Box test

struct LjungBoxTest <: HypothesisTest
    y::Vector{Float64}  # time series vector
    n::Int              # number of observations
    lag::Int            # number of lags in test statistic
    dof::Int            # degrees-of-freedom correction (if y is residual)
    Q::Float64          # test statistic
end

"""
    LjungBoxTest(y, lag, dof=0)

Compute the Ljung-Box `Q` statistic to test the null hypothesis of
independence in a time series `y`.

`lag` specifies the number of lags used in the construction of `Q`. When
testing the residuals of an estimated model, `dof` has to be set to the number
of estimated parameters. E.g., when testing the residuals of an ARIMA(p,0,q)
model, set `dof=p+q`.

# External links

  * [Ljung-Box test on Wikipedia](https://en.wikipedia.org/wiki/Ljung–Box_test)
"""
function LjungBoxTest(y::AbstractVector{T}, lag::Int, dof::Int=0) where T<:Real
    if dof>=lag
        throw(ArgumentError("Number of lags has to be larger than degrees of" *
        " freedom correction"))
    end
    n = size(y,1)
    Q = n*(n+2)* first(sum(k -> autocor(y,k:k).^2/(n-k), 1:lag))
    LjungBoxTest(y,n,lag,dof,Q)
end

testname(::LjungBoxTest) = "Ljung-Box autocorrelation test"
population_param_of_interest(x::LjungBoxTest) = ("autocorrelations up to lag k",
    "all zero", NaN)
default_tail(test::LjungBoxTest) = :right

function show_params(io::IO, x::LjungBoxTest, ident)
    println(io, ident, "number of observations:         ", x.n)
    println(io, ident, "number of lags:                 ", x.lag)
    println(io, ident, "degrees of freedom correction:  ", x.dof)
    println(io, ident, "Q statistic:                    ", x.Q)
end

pvalue(x::LjungBoxTest) = pvalue(Chisq(x.lag-x.dof), x.Q; tail=:right)
