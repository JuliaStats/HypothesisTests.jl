# white.jl
# White and Breusch-Pagan tests for heteroskedasticity
#
# Copyright (C) 2020   Paul Söderlind
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

export WhiteTest, BreuschPaganTest

struct WhiteTest <: HypothesisTest
    dof::Int                 #degrees of freedom in test
    LM::Float64              #test statistic, distributed as Chisq(dof)
    TestDescription::String  #short description of the test
end

"""
    WhiteTest(X, e, TestType)

Compute White's (or Breusch-Pagan's) test for heteroskedasticity.

`X` is a matrix of regressors and `e` is the vector of residuals from the original model.
`TestType` is either `:Linear` for the Breusch-Pagan/Koenker test, `:LinearAndSquares` for
White's test with linear and squared terms only (no cross-products), or `:White` (the default)
for the full White's test (linear, squared and cross-product terms). `X` should include a
constant. In some applications, `X` is a subset of the regressors in the original model,
or just the fitted values. This saves degrees of freedom and may give a more powerful test.
The test statistic is T*R2 where R2 is from the regression of `e^2` on the terms
mentioned above. Under the null hypothesis it is distributed as `Chisquare(dof)`
where `dof` is the number of independent terms (not counting the constant), so the null
is rejected when the test statistic is large enough.

Implements: [`pvalue`](@ref)

# References
  * H. White, (1980): A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity, Econometrica, 48, 817-838.
  * T.S. Breusch & A.R. Pagan (1979), A simple test for heteroscedasticity and random coefficient variation, Econometrica, 47, 1287-1294
  * R. Koenker (1981), A note on studentizing a test for heteroscedasticity, Journal of Econometrics, 17, 107-112

# External links
  * [White's test on Wikipedia](https://en.wikipedia.org/wiki/White_test)
"""
function WhiteTest(X::AbstractVecOrMat{<:Real}, e::AbstractVector{<:Real}, TestType=:White)
    (n,K) = (size(X,1),size(X,2))    #nobs,nvars

    n == length(e) || throw(DimensionMismatch("inputs must have the same length"))

    if TestType == :Linear
        z  = X
        TestDescription = "Breusch-Pagan's test for heteroskedasticity (linear)"
    elseif TestType == :LinearAndSquares
        z = [X X.^2]
        TestDescription = "White's test for heteroskedasticity (linear and squares)"
    else               #linear, squares and cross-products
        z = fill(NaN,n,round(Int,K*(K+1)/2))   #loop is easier and quicker than index tricks
        vv = 1
        for i = 1:K, j = i:K
            z[:,vv] = X[:,i].*X[:,j]             #eg. x1*x1, x1*x2, x2*x2
            vv      = vv + 1
        end
        TestDescription = "White's test for heteroskedasticity (linear, squares, cross products)"
    end

    dof  = rank(z) - 1                         #number of independent regressors in z
    e2  = e.^2
    b   = z\e2
    res = e2 - z*b
    R2  = 1 - var(res)/var(e2)
    LM  = n*R2                                 #n*R2 in original papers, could also use n*R2/(1-R2)
    return WhiteTest(dof, LM, TestDescription)
end

"""
    BreuschPaganTest(X, e)

Compute Breusch-Pagan's test for heteroskedasticity.

`X` is the matrix of regressors from the original model and `e` the vector of residuals.
See `WhiteTest` for further details.
"""
BreuschPaganTest(X, e) = WhiteTest(X, e, :Linear)


testname(t::WhiteTest)                     = t.TestDescription
population_param_of_interest(t::WhiteTest) = ("T*R2", 0, t.LM)
default_tail(test::WhiteTest)              = :right

function show_params(io::IO, t::WhiteTest, ident)
    println(io, ident, "T*R^2 statistic:        ", t.LM)
    println(io, ident, "degrees of freedom:     ", t.dof)
end

StatsBase.dof(t::WhiteTest) = t.dof
pvalue(t::WhiteTest) = pvalue(Chisq(t.dof), t.LM, tail=:right)
