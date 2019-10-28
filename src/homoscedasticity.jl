# homoscedasticity.jl
# Homoscedasticity test
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

export HomoscedasticityTest

struct HomoscedasticityTest <: HypothesisTest
    n::Int        # number of observations
    h::Int        # separator
    H::Float64    # test statistic
end

"""
    HomoscedasticityTest(y::AbstractVector, h::Int)

Compute the test statistic of an homoscedasticity test: the null hypothesis is that a real-valued
vector `y` has constant variance.

# References

  * James Durbin and Siem Jan Koopman, "Time Series Analysis by State Space Methods", 2nd Ed.,
    2012, Oxford Statistical Science Series.

# External links

  * [Homoscedasticity on Wikipedia](https://en.wikipedia.org/wiki/Homoscedasticity)
"""
function HomoscedasticityTest(y::AbstractVector{T}, h::Int) where T <: Real
    n = length(y)
    @assert h < n/2
    H = sum(y[t]^2 for t = n-h+1:n)/sum(y[t]^2 for t = 1:h)
    return HomoscedasticityTest(n, h, H)
end

testname(::HomoscedasticityTest) = "Homoscedasticity test"
population_param_of_interest(x::HomoscedasticityTest) =
    ("variance ratio", "1.0", "$(x.H)")
default_tail(test::HomoscedasticityTest) = :both

function show_params(io::IO, x::HomoscedasticityTest, ident)
    println(io, ident, "number of observations:         ", x.n)
    println(io, ident, "H statistic:                    ", x.H)
end

function pvalue(x::HomoscedasticityTest; tail = :both)
    dist = FDist(x.h, x.h)
    if tail == :both
        Δ = abs(cdf(dist, x.H) - 0.5)
        H1 = quantile(dist, 0.5 - Δ)
        H2 = quantile(dist, 0.5 + Δ)
        return cdf(dist, H1) + (1 - cdf(dist, H2))
    elseif tail == :right
        return 1 - cdf(dist, x.H)
    elseif tail == :left
        return cdf(dist, x.H)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end