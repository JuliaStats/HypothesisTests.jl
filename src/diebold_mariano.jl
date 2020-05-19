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
    df::Int        # degrees of freedom
    t::Float64     # test statistic
end

"""
    DieboldMarianoTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; loss=abs2, h=1)

Perform the modified Diebold-Mariano test proposed by Harvey, Leybourne and Newbold of the null 
hypothesis that the two methods have the same forecast accuracy. `loss` is the loss function described
in Diebold, F.X. and Mariano, R.S. (1995) Comparing predictive accuracy. Journal of Business and 
Economic Statistics, 13, 253-263. and `h` is the number of steps ahead of the forecast.

Implements: [`pvalue`](@ref)

# References

  * Diebold, F.X. and Mariano, R.S. (1995) Comparing predictive accuracy. 
    Journal of Business and Economic Statistics, 13, 253-263.

  * Harvey, D., Leybourne, S., & Newbold, P. (1997). Testing the equality of prediction 
    mean squared errors. International Journal of forecasting, 13(2), 281-291.
  
"""
function DieboldMarianoTest(e1::AbstractVector{<:Real}, e2::AbstractVector{<:Real}; 
                            loss::Function=abs2, h::Integer=1)

    @assert length(e1) == length(e2)
    n = length(e1)
    # Calculate the loss diferential series based on the loss function g
    d = loss.(e1) .- loss.(e2)
    dm_cov = autocov(d, collect(0:h - 1))
    dm_var = sum([dm_cov[1]; 2 * dm_cov[2:end]])/n
    # Statistic from the original Diebold-Mariano test 
    statistic_dm = mean(d)/sqrt(dm_var)
    k = sqrt((1 + (1 - 2*h + (h/n)*(h - 1))/n))
    # Statistic from the modified Diebold-Mariano test proposed by Harvey, Leybourne and Newbold
    statistic_hln = statistic_dm * k
    return DieboldMarianoTest(n, n - 1, statistic_hln)
end

testname(::DieboldMarianoTest) = "Diebold Mariano test"
population_param_of_interest(x::DieboldMarianoTest) = ("mean", 0.0, x.t)
default_tail(test::DieboldMarianoTest) = :both

function show_params(io::IO, x::DieboldMarianoTest, ident)
    println(io, ident, "number of observations: $(x.n)")
    println(io, ident, "DM statistic:           $(x.t)")
    println(io, ident, "degrees of freedom:     $(x.df)")
end
