# f.jl
# F-tests
#
# Copyright (C) 2019   Raphael Saavedra
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

export VarianceFTest

struct VarianceFTest <: HypothesisTest
    n_x::Int   # size of sample 1
    n_y::Int   # size of sample 2
    df_x::Int  # degrees of freedom 1
    df_y::Int  # degrees of freedom 2
    F::Real    # test statistic
    F0::Real   # variance ratio under h_0
end

"""
    VarianceFTest(y1::AbstractVector{<: Real}, y2::AbstractVector{<: Real})

Compute the test statistic of an F test: the null hypothesis is that two real-valued vectors
`y1` and `y2` have equal variances.

# References

  * George E. P. Box, "Non-Normality and Tests on Variances", Biometrika 40 (3/4): 318–335, 1953.

# External links

  * [F-test of equality of variances on Wikipedia](https://en.wikipedia.org/wiki/F-test_of_equality_of_variances)
"""
function VarianceFTest(y1::AbstractVector{<: Real}, y2::AbstractVector{<: Real}; F0::Real = 0)
    n1, n2 = length(y1), length(y2)
    F = var(y1) / var(y2)
    return VarianceFTest(n1, n2, n1-1, n2-1, F, F0)
end

"""
    VarianceFTest(y::AbstractVector{<: Real}, h::Integer)

Compute the test statistic of an homoscedasticity F test: the null hypothesis is that a real-valued
vector `y` has constant variance. The test is conducted over the ratio of the sample variances of the
first `h` observations and the last `h` observations of the series.

# References

  * James Durbin and Siem Jan Koopman, "Time Series Analysis by State Space Methods", 2nd Ed.,
    2012, Oxford Statistical Science Series.

# External links

  * [Homoscedasticity on Wikipedia](https://en.wikipedia.org/wiki/Homoscedasticity)
"""
function VarianceFTest(y::AbstractVector{<: Real}, h::Integer)
    n = length(y)
    h > n/2 && throw(ArgumentError("The number of observations considered in each end must not be more than half the total number of observations"))
    y1 = y[n-h+1:n]
    y2 = y[1:h]
    return VarianceFTest(y1, y2)
end

testname(::VarianceFTest) = "Variance F-test"
population_param_of_interest(x::VarianceFTest) = ("variance ratio", x.F0, x.F)
default_tail(test::VarianceFTest) = :both

function show_params(io::IO, x::VarianceFTest, ident)
    println(io, ident, "number of observations: [$(x.n_x), $(x.n_y)]")
    println(io, ident, "F statistic:            $(x.F)")
    println(io, ident, "degrees of freedom:     [$(x.df_x), $(x.df_y)]")
end

function pvalue(x::VarianceFTest; tail = :both)
    dist = FDist(x.df_x, x.df_y)
    if tail == :both
        Δ = abs(cdf(dist, x.F) - 0.5)
        F1 = quantile(dist, 0.5 - Δ)
        F2 = quantile(dist, 0.5 + Δ)
        return cdf(dist, F1) + (1 - cdf(dist, F2))
    elseif tail == :right
        return 1 - cdf(dist, x.F)
    elseif tail == :left
        return cdf(dist, x.F)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end