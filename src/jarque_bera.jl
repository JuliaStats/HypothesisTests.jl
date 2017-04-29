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

immutable JarqueBeraTest <: HypothesisTest
    n       ::Int           # number of observations
    JB      ::Float64       # test statistic
    m3      ::Float64       # skewness
    m4      ::Float64       # excess kurtosis
end

"""
    JarqueBeraTest(y)

Compute the Jarque-Bera statistic to test the null hypothesis that a series `y` comes from a normal distribution.

Note that the p-values might be distorted in small samples.

External links

* [Jarque-Bera test on Wikipedia](https://en.wikipedia.org/wiki/Jarque–Bera_test)
"""

function JarqueBeraTest{T<:Real}(y::AbstractVector{T})
    n = length(y)
    m1 = sum(y)/n
    m2 = sum((y - m1).^2)/n
    m3 = sum((y - m1).^3)/n
    m4 = sum((y - m1).^4)/n
    b1 = (m3/m2^(3/2))^2
    b2 = (m4/m2^2)

    stat = n * b1/6 + n*(b2 - 3)^2/24
    JarqueBeraTest(n, stat, m3, m4)
end

testname(::JarqueBeraTest) = "Jarque-Bera normality test"
population_param_of_interest(x::JarqueBeraTest) =
    ("skewness and excess kurtosis", "both zero", "$(x.m3) and $(x.m4-3)")
default_tail(test::JarqueBeraTest) = :right

function show_params(io::IO, x::JarqueBeraTest, ident)
    println(io, ident, "number of observations:         ", x.n)
    println(io, ident, "JB statistic:                   ", x.JB)
end

pvalue(x::JarqueBeraTest) = pvalue(Chisq(2), x.JB; tail=:right)
