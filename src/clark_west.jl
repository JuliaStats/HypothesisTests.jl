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

struct ClarkWestTest <: TTest
    n::Integer     # number of observations
    df::Integer    # degrees of freedom
    t::Real        # test statistic
end

"""
    ClarkWestTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Perform the Clark-West test ...

Implements: [`pvalue`](@ref)

# References

 * Clark, T. E., West, K. D. 2006, Using out-of-sample mean squared prediction errors to test 
   the martingale difference hypothesis. Journal of Econometrics, 135(1): 155–186.

 * Clark, T. E., West, K. D. 2007, Approximately normal tests for equal predictive accuracy 
   in nested models. Journal of Econometrics, 138(1): 291–311. 
"""
function ClarkWestTest(e1::AbstractVector{<:Real}, e2::AbstractVector{<:Real}; 
                        h::Integer=1)

    @assert length(e1) == length(e2)
    n = length(e1)
    
    
    
    return ClarkWestTest(n, n - 1, statistic_cw)
end

testname(::ClarkWestTest) = "Clark West test"
population_param_of_interest(x::ClarkWestTest) = ("mean", 0.0, x.t)
default_tail(test::ClarkWestTest) = :both

function show_params(io::IO, x::ClarkWestTest, ident)
    println(io, ident, "number of observations: $(x.n)")
    println(io, ident, "CW statistic:           $(x.t)")
    println(io, ident, "degrees of freedom:     $(x.df)")
end