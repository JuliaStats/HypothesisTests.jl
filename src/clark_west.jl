# clark_west.jl
# Clark-West
#
# Copyright (C) 2020   Guilherme Bodin, Paul Söderlind
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

export ClarkWestTest

struct ClarkWestTest <: ZTest
    n::Int       # number of observations
    xbar::Real   # estimated mean
    stderr::Real # population standard error
    z::Real      # t-statistic
    μ0::Real     # mean under h_0
end

"""
    ClarkWestTest(e1::AbstractVector{<:Real}, e2::AbstractVector{<:Real}, lookahead::Integer=1)

Perform the Clark-West test of equal performance of two nested prediction models, in terms of the
out-of-sample mean squared prediction errors.

`e1` is a vector of forecasts from the smaller (nested) model, `e2` is a vector of forecast
errors from the larger model, and `lookahead` is the number of steps ahead of the forecast.
Typically, the null hypothesis is that the two models perform equally well (a two-sided test),
but sometimes we test whether the larger model performs better, which is indicated by a
positive test statistic, for instance, above 1.645 for the 5% significance level (right tail test).

Implements: [`pvalue`](@ref)
# References
 * Clark, T. E., West, K. D. 2006, Using out-of-sample mean squared prediction errors to test
   the martingale difference hypothesis. Journal of Econometrics, 135(1): 155–186.
 * Clark, T. E., West, K. D. 2007, Approximately normal tests for equal predictive accuracy
   in nested models. Journal of Econometrics, 138(1): 291–311.

"""
function ClarkWestTest(e1::AbstractVector{<:Real}, e2::AbstractVector{<:Real},
                       lookahead::Integer=1)
    length(e1) == length(e2) || throw(DimensionMismatch("inputs must have the same length"))
    n            = length(e1)
    d            = 2*e1.*(e1 - e2)
    cw_cov       = autocov(d, 0:lookahead-1)
    cw_var       = (cw_cov[1] + 2*sum(@view(cw_cov[2:end])))/n
    xbar         = mean(d)
    stderr       = sqrt(cw_var)
    statistic_cw = xbar/stderr
    return ClarkWestTest(n, xbar, stderr, statistic_cw, 0)
end

testname(::ClarkWestTest) = "Clark West test"
population_param_of_interest(x::ClarkWestTest) = ("Mean", x.μ0, x.xbar)
default_tail(::ClarkWestTest) = :both

function show_params(io::IO, x::ClarkWestTest, ident)
    println(io, ident, "number of observations:    $(x.n)")
    println(io, ident, "CW statistic:              $(x.z)")
    println(io, ident, "population standard error: $(x.stderr)")
end
