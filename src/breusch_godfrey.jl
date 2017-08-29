# breusch_godfrey.jl
# Breusch-Godfrey test for autocorrelation
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

export BreuschGodfreyTest

struct BreuschGodfreyTest <: HypothesisTest
    n::Int              # number of observations
    lag::Int            # number of lags in test statistic
    BG::Float64         # test statistic
end

"""
    BreuschGodfreyTest(X, e, lag, start0 = true)

Compute the Breusch-Godfrey test for serial correlation in the residuals
of a regression model.

`X` is the matrix of regressors from the original model and `e` the vector of residuals.
`lag` determines the number of lagged residuals included in the auxiliary regression.
Set `start0` to specify how the starting values for the lagged residuals are handled.
`start0 = true` (default) sets them to zero (as in Godfrey, 1978); `start0 = false`
uses the first `lag` residuals as starting values, i.e. shortening the sample by `lag`.

# External links

  * [Breusch-Godfrey test on Wikipedia](https://en.wikipedia.org/wiki/Breusch–Godfrey_test)
"""
function BreuschGodfreyTest(xmat::AbstractArray{T}, e::AbstractVector{T},
                            lag::Int, start0::Bool=true) where T<:Real
    n = size(e,1)
    elag = zeros(Float64,n,lag)
    for ii = 1:lag  # construct lagged residuals
        elag[ii+1:end,ii] = e[1:end-ii]
    end

    offset = start0 ? 0 : lag

    regmat = [xmat[offset+1:end,:] elag[offset+1:end,:]]
    regcoeff = regmat\e[offset+1:end]
    resid = e[offset+1:end] - regmat*regcoeff

    rsq = 1 - dot(resid,resid)/dot(e[offset+1:end],e[offset+1:end]) # uncentered R^2
    BG = (n-offset)*rsq
    BreuschGodfreyTest(n-offset,lag,BG)
end

testname(::BreuschGodfreyTest) = "Breusch-Godfrey autocorrelation test"
population_param_of_interest(x::BreuschGodfreyTest) =
    ("coefficients on lagged residuals up to lag p", "all zero", NaN)
default_tail(test::BreuschGodfreyTest) = :right

function show_params(io::IO, x::BreuschGodfreyTest, ident)
    println(io, ident, "number of observations:         ", x.n)
    println(io, ident, "number of lags:                 ", x.lag)
    println(io, ident, "T*R^2 statistic:                ", x.BG)
end

pvalue(x::BreuschGodfreyTest) = pvalue(Chisq(x.lag), x.BG; tail=:right)
