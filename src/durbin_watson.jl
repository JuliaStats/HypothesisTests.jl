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

By default, the choice of approach to compute p-values depends on the sample size
(`p_compute=:ndep`). For small samples (n<100), Pan's algorithm (Farebrother, 1980) is
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

  * [Durbin-Watson test on Wikipedia](https://en.wikipedia.org/wiki/Durbin–Watson_statistic)
"""
function DurbinWatsonTest(xmat::AbstractArray{T}, e::AbstractArray{T};
    p_compute::Symbol = :ndep) where T<:Real

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

"""
    pan_algorithm(a::AbstractArray, x::Float64, m::Int, n::Int)

Compute exact p-values for the Durbin-Watson statistic using Pan's algorithm (Farebrother,
1980).

`a` is the vector of non-zero Eigenvalues of ``(I-(X(X'X)^{-1}X'))A`` (see Durbin and
Watson, 1971, p. 2), `x` is the value of the Durbin-Watson statistic, `m` the number of
elements in `a`, and `n` the number of approximation terms (see Farebrother, 1980, eq. 5).

# References

  * J. Durbin and G. S. Watson, 1971, "Testing for Serial Correlation in Least Squares
  Regression: III", Biometrika, Vol. 58, No. 1, pp. 1-19,
  [http://www.jstor.org/stable/2334313](http://www.jstor.org/stable/2334313).
  * R. W. Farebrother, 1980, "Algorithm AS 153: Pan's Procedure for the Tail Probabilities
  of the Durbin-Watson Statistic", Journal of the Royal Statistical Society, Series C
  (Applied Statistics), Vol. 29, No. 2, pp. 224-227,
  [http://www.jstor.org/stable/2986316](http://www.jstor.org/stable/2986316).

"""
function pan_algorithm(a::AbstractArray, x::Float64, m::Int, n::Int)

    ν = findfirst(ai -> ai >= x, a)
    if ν == 0 || ν === nothing
        return 1.0
    elseif ν == 1
        return 0.0
    else
        k = 1
        ν = ν - 1
        h = m - ν

        if ν <= h
            d  = 2; h  = ν; k  = - k; j1 = 0; j2 = 2; j3 = 3; j4 = 1
        else
            d  = - 2; ν  = ν + 2; j1 = m - 2; j2 = m - 1; j3 = m + 1; j4 = m
        end

        pin = pi / (2n)
        sum = (k + 1) / 2
        sgn = k / n
        n2  = 2n - 1

        # first integral
        for  f1 = h - 2 * floor(Int,h/2) : -1 : 0
            for f2 = j2:d:ν
                sum1 = a[j4]
                if f2 == 0
                    prod = x
                else
                    prod = a[f2]
                end
                u = 0.5 * (sum1 + prod)
                v = 0.5 * (sum1 - prod)
                sum1 = 0.0
                for  i = n2:-2:1
                    y = u - v * cos(i * pin)
                    num = y - x
                    prod = 1.0
                    for k in [(1:j1)' (j3:m)']
                        prod *= num / (y - a[k])
                    end
                    sum1 += sqrt(abs(prod))
                end
                sgn = -sgn
                sum += sgn * sum1
                j1 += d
                j3 += d
                j4 += d
            end

            # second integral
            if d == 2
                j3 = j2
            else
                j1 = j2
            end
            j2 = 0
            ν  = 0
        end
        return sum
    end

end

function pvalue(x::DurbinWatsonTest; tail=:both)

    exact_problem_flag = 0
    if (x.p_compute == :ndep && x.n <= 100) || x.p_compute == :exact
        # p-vales based on Pan's algorithm (see Farebrother, 1980)

        # the following setup is, e.g, described in Durbin and Watson (1971)
        A = diagm(-1 => -ones(x.n - 1), 0 => 2 * ones(x.n), 1 => -ones(x.n - 1))
        A[1, 1] = 1
        A[x.n, x.n] = 1
        EV_temp = sort(real(eigvals(A - x.xmat / (x.xmat' * x.xmat) * (x.xmat'*A))))
        EV = EV_temp[EV_temp .> 1e-10]
        p_temp = pan_algorithm(EV, x.DW, length(EV), 15)

        if p_temp < 0.0 || p_temp > 1.0
            # println(p_temp)
            @warn("Exact p-values outside [0,1]. Approximate p-values reported instead.")
            exact_problem_flag = 1
        end
    end

    if exact_problem_flag == 1 || (x.p_compute == :ndep && x.n > 100
        ) || x.p_compute == :approx
        # p-values based on normal approximation (see Durbin and Watson, 1950)

        # the following derivations follow Durbin and Watson (1951, p. 164)
        X = x.xmat
        k = size(X, 2)
        inv_XX = inv(X' * X)

        AX = zeros(x.n, k)
        AX[[1, x.n], :] = X[[1, x.n], :] - X[[2, x.n - 1], :]
        for i = 2:(x.n - 1)
            AX[i, :] = - X[i - 1, :] + 2 * X[i, :] - X[i + 1, :]
        end

        temp_mat = X' * AX * inv_XX
        P = 2 * (x.n - 1) - tr(temp_mat) # first term: tr(A)
        Q = 2 * (3 * x.n - 4) - 2 * tr(AX' * AX * inv_XX) + tr(temp_mat^2)
        # first term: tr(A^2)

        dw_mean = P / (x.n - k)
        dw_var = 2 / ((x.n - k) * (x.n - k + 2)) * (Q - P * dw_mean)

        p_temp = cdf(Normal(dw_mean, sqrt(dw_var)), x.DW)
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
