# diebold_mariano.jl
# Diebold-Mariano
#
# Copyright (C) 2020   Guilherme Bodin
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

export DieboldMarianoTest

struct DieboldMarianoTest <: TTest
    n::Int         # number of observations
    xbar::Real     # estimated mean
    df::Int        # degrees of freedom
    t::Real        # test statistic
    stderr::Real   # empirical standard error
    μ0::Real       # mean under h_0
end

"""
    DieboldMarianoTest(e1::AbstractVector{<:Real}, e2::AbstractVector{<:Real}; loss=abs2, lookahead=1)

Perform the modified Diebold-Mariano test proposed by Harvey, Leybourne and Newbold of the null 
hypothesis that the two methods have the same forecast accuracy. `loss` is the loss function described
in Diebold, F.X. and Mariano, R.S. (1995) Comparing predictive accuracy. Journal of Business and 
Economic Statistics, 13, 253-263. and `lookahead` is the number of steps ahead of the forecast.

# References

  * Diebold, F.X. and Mariano, R.S. (1995) Comparing predictive accuracy. 
    Journal of Business and Economic Statistics, 13, 253-263.

  * Harvey, D., Leybourne, S., & Newbold, P. (1997). Testing the equality of prediction 
    mean squared errors. International Journal of forecasting, 13(2), 281-291.
  
"""
function DieboldMarianoTest(e1::AbstractVector{<:Real}, e2::AbstractVector{<:Real}; 
                            loss::Function=abs2, lookahead::Integer=1)
    length(e1) == length(e2) || throw(DimensionMismatch("inputs must have the same length"))
    n = length(e1)
    # Calculate the loss diferential series based on the loss function g
    d = loss.(e1) .- loss.(e2)
    dm_cov = autocov(d, collect(0:lookahead-1))
    dm_var = (dm_cov[1] + 2 * sum(dm_cov[2:end]))/n
    # Statistic from the original Diebold-Mariano test 
    xbar_dm = mean(d)
    stderr_dm = sqrt(dm_var)
    statistic_dm = xbar_dm/stderr_dm
    k = sqrt((1 + (1 - 2*lookahead + (lookahead/n)*(lookahead - 1))/n))
    # Statistic from the modified Diebold-Mariano test proposed by Harvey, Leybourne and Newbold
    statistic_hln = statistic_dm * k
    xbar = xbar_dm * k
    stderr = stderr_dm
    return DieboldMarianoTest(n, xbar, n - 1, statistic_hln, stderr, 0.0)
end

testname(::DieboldMarianoTest) = "Diebold-Mariano test"
population_param_of_interest(x::DieboldMarianoTest) = ("Mean", x.μ0, x.xbar)
default_tail(test::DieboldMarianoTest) = :both

function show_params(io::IO, x::DieboldMarianoTest, ident)
    println(io, ident, "number of observations: $(x.n)")
    println(io, ident, "DM statistic:           $(x.t)")
    println(io, ident, "degrees of freedom:     $(x.df)")
end
