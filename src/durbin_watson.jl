# durbin_watson.jl
# Durbin-Watson test for autocorrelation
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

export DurbinWatsonTest

struct DurbinWatsonTest <: HypothesisTest
    xmat::Array{Float64}  # regressor matrix
    n::Int                # number of observations
    DW::Float64           # test statistic
    p_compute::Symbol     # determines how p-values are computed
end

"""
    DurbinWatsonTest(X::AbstractArray, e::AbstractVector; p_compute::Symbol = :ndep)

Compute the Durbin-Watson test for serial correlation in the residuals of a regression model.

`X` is the matrix of regressors from the original regression model and `e` the vector of
residuals. Note that the Durbin-Watson test is not valid if `X` includes a lagged dependent
variable. The test statistic is computed as
```math
DW = \\frac{\\sum_{t=2}^n (e_t - e_{t-1})^2}{\\sum_{t=1}^n e_t^2}
```
where `n` is the number of observations.

By default, the choice of approach to compute p-values depends on the sample size (`p_compute=:ndep`). For small samples (n<100), Pan's algorithm (Farebrother, 1980) is
employed. For larger samples, a normal approximation is used (Durbin and Watson, 1950). To
always use Pan's algorithm, set `p_compute=:exact`. `p_compute=:approx` will always use the
normal approximation.

Default is a two-sided p-value for the alternative hypothesis of positive or negative
serial correlation. One-sided p-values can be requested by calling
`pvalue(x::DurbinWatsonTest; tail=)` with the options `:left` (negative serial correlation)
and `:right` (positive serial correlation).

# References

  * J. Durbin and G. S. Watson, 1951, "Testing for Serial Correlation in Least Squares
  Regression: II", Biometrika, Vol. 38, No. 1/2, pp. 159-177,
  [http://www.jstor.org/stable/2332325](http://www.jstor.org/stable/2332325).
  * J. Durbin and G. S. Watson, 1950, "Testing for Serial Correlation in Least Squares
  Regression: I", Biometrika, Vol. 37, No. 3/4, pp. 409-428,
  [http://www.jstor.org/stable/2332391](http://www.jstor.org/stable/2332391).
  * R. W. Farebrother, 1980, "Algorithm AS 153: Pan's Procedure for the Tail Probabilities
  of the Durbin-Watson Statistic", Journal of the Royal Statistical Society, Series C
  (Applied Statistics), Vol. 29, No. 2, pp. 224-227,
  [http://www.jstor.org/stable/2986316](http://www.jstor.org/stable/2986316).

# External links

  * [Durbin-Watson test on Wikipedia:
  https://en.wikipedia.org/wiki/Durbin–Watson_statistic
  ](https://en.wikipedia.org/wiki/Durbin–Watson_statistic)
"""
function DurbinWatsonTest{T<:Real}(xmat::AbstractArray{T}, e::AbstractVector{T};
    p_compute::Symbol = :approx) # approx as ndep not implemented completely

    n = length(e)
    DW = sum(diff(e) .^2) / sum(e .^2)

    DurbinWatsonTest(xmat, n, DW, p_compute)
end

testname(::DurbinWatsonTest) = "Durbin-Watson autocorrelation test"
population_param_of_interest(x::DurbinWatsonTest) =
    ("sample autocorrelation parameter", "0", 1 - x.DW / 2)
default_tail(test::DurbinWatsonTest) = :both

function show_params(io::IO, x::DurbinWatsonTest, ident)
    println(io, ident, "number of observations:     ", x.n)
    println(io, ident, "DW statistic:               ", x.DW)
end

function pvalue(x::DurbinWatsonTest; tail=:both)

    if (x.p_compute == :ndep && x.n > 100) || x.p_compute == :approx
        # p-values based on normal approximation (see Durbin and Watson, 1950)

        # the following derivations follow Durbin and Watson (1951, p. 164)
        X = x.xmat
        k = size(X, 2)
        inv_XX = (X' * X) \ eye(k)

        AX = zeros(x.n, k)
        AX[[1, x.n], :] = X[[1, x.n], :] - X[[2, x.n - 1], :]
        for i = 2:(x.n - 1)
            AX[i, :] = - X[i - 1, :] + 2 * X[i, :] - X[i + 1, :]
        end

        temp_mat = X' * AX * inv_XX
        P = 2 * (x.n - 1) - trace(temp_mat) # first term: trace(A)
        Q = 2 * (3 * x.n - 4) - 2 * trace(AX' * AX * inv_XX) + trace(temp_mat^2)
        # first term: trace(A^2)

        dw_mean = P / (x.n - k)
        dw_var = 2 / ((x.n - k) * (x.n - k + 2)) * (Q - P * dw_mean)

        p_temp = cdf(Normal(dw_mean, sqrt(dw_var)), x.DW)

    elseif (x.p_compute == :ndep && x.n <= 100) || x.p_compute == :exact
        # p-vales based on Pan's algorithm (see Farebrother, 1980)
        throw(ArgumentError("not implemented yet"))
    end

    if tail == :both
        2 * min(p_temp, 1 - p_temp)
    elseif tail == :right
        p_temp
    elseif tail == :left
        1 - p_temp
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end

end
