# jarque_bera.jl
# Jarque-Bera goodness-of-fit test
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

export JarqueBeraTest

struct JarqueBeraTest <: HypothesisTest
    n::Int         # number of observations
    JB::Float64    # test statistic
    skew::Float64  # skewness
    kurt::Float64  # excess kurtosis
end

"""
    JarqueBeraTest(y::AbstractVector; adjusted::Bool=false)

When `adjusted` is `false`, compute the Jarque-Bera statistic to test the null hypothesis
that a real-valued vector `y` is normally distributed.

Note that the approximation by the Chi-squared distribution does not work well and the
speed of convergence is slow. In small samples, the test tends to be over-sized for
nominal levels up to about 3% and under-sized for larger nominal levels (Mantalos, 2010).

When `adjusted` is `true`, compute the Adjusted Lagrangian Multiplier statistic to test the
null hypothesis that a real-valued vector `y` is normally distributed.

Note that the use of Adjusted Lagrangian Multiplier is preferred over Jarque-Bera for small
and medium sample sizes and it is a modification to the Jarque-Bera test (Urzua, 1996).

# References

  * Panagiotis Mantalos, 2011, "The three different measures of the sample skewness and
    kurtosis and the effects to the Jarque-Bera test for normality", International Journal
    of Computational Economics and Econometrics, Vol. 2, No. 1,
    [link](http://dx.doi.org/10.1504/IJCEE.2011.040576).

  * Carlos M. Urzúa, "On the correct use of omnibus tests for normality", Economics Letters,
    Volume 53, Issue 3,
    [link](https://doi.org/10.1016/S0165-1765(96)00923-8).

# External links

  * [Jarque-Bera test on Wikipedia](https://en.wikipedia.org/wiki/Jarque–Bera_test)
"""
function JarqueBeraTest(y::AbstractVector{T}; adjusted::Bool=false) where T<:Real
    n = length(y)
    M = Base.promote_op(/, T, typeof(n))
    m1r = m2r = m3r = m4r = zero(M)
    @inbounds for yi in y # compute raw moments
        m1r += yi / n
        m2r += yi^2 / n
        m3r += yi^3 / n
        m4r += yi^4 / n
    end
    # compute central moments (http://mathworld.wolfram.com/CentralMoment.html)
    m2 = -m1r^2 + m2r
    m3 = 2 * m1r^3 - 3 * m1r * m2r + m3r
    m4 = -3 * m1r^4 + 6 * m1r^2 * m2r - 4 * m1r * m3r + m4r

    skew = m3 / m2^(3/2)
    kurt = m4 / m2^2

    if adjusted == false
        stat = n * skew^2 / 6 + n * (kurt - 3)^2 / 24
    else
        meankurt = 3 * (n-1) / (n+1)
        varskew = 6 * (n-2) / ((n+1) * (n+3))
        varkurt = 24 * n * (n-2) * (n-3) / ((n+1)^2 * (n+3) * (n+5))

        stat = skew^2 / varskew + (kurt - meankurt)^2 / varkurt
    end

    JarqueBeraTest(n, stat, skew, kurt)
end

testname(::JarqueBeraTest) = "Jarque-Bera normality test"
population_param_of_interest(x::JarqueBeraTest) =
    ("skewness and kurtosis", "0 and 3", "$(x.skew) and $(x.kurt)")
default_tail(test::JarqueBeraTest) = :right

function show_params(io::IO, x::JarqueBeraTest, ident)
    println(io, ident, "number of observations:         ", x.n)
    println(io, ident, "JB statistic:                   ", x.JB)
end

pvalue(x::JarqueBeraTest) = pvalue(Chisq(2), x.JB; tail=:right)
