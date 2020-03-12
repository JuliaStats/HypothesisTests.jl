# augmented_dickey_fuller.jl
# (Augmented) Dickey-Fuller unit root test
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

export ADFTest

struct ADFTest <: HypothesisTest
    n               ::Int               # number of observations
    deterministic   ::Symbol            # deterministic terms included in regression
    lag             ::Int               # number of lags in test statistic
    stat            ::Float64           # test statistic
    coef            ::Float64           # coefficient on lagged (non-differenced) variable
    cv              ::Vector{Float64}   # critical values
end

"""
    ADFTest(y::AbstractVector{T}, deterministic::Symbol, lag::Int) where T<:Real

Compute the augmented Dickey-Fuller unit root test.

`y` is the time series to be tested, `deterministic` determines the deterministic terms
(options: `:none`, `:constant`, `:trend`, `:squared_trend`) and `lag` the number of lagged
first-differences included in the test regression, respectively.

Critical values and asymptotic p-values are computed based on response surface regressions
following MacKinnon (2010) and MacKinnon (1994), respectively. These may differ slightly
from those reported in other regression packages as different algorithms might be used.

# References

  * James G. MacKinnon, 2010, "Critical values for cointegration tests,"
    QED Working Paper No. 1227, 2010, [link](http://ideas.repec.org/p/qed/wpaper/1227.html).
  * James G. MacKinnon, 1994, "Approximate Asymptotic Distribution Functions for
    Unit-Root and Cointegration Tests", Journal of Business & Economic Statistics,
    Vol. 12, No. 2, pp. 167-176, [link](http://www.jstor.org/stable/1391481).

# External links

  * [Augmented Dickey-Fuller test on Wikipedia](https://en.wikipedia.org/wiki/Augmented_Dickey–Fuller_test)
"""
function ADFTest(y::AbstractVector{T}, deterministic::Symbol, lag::Int) where T<:Real

    nobs = length(y)
    Δy = diff(y)

    if deterministic == :none
        Δyt   = Δy[lag+1:end]
        xt    = y[lag+1:end-1,1:1]
        dfexo = 0
    elseif deterministic == :constant
        Δyt   = Δy[lag+1:end]
        xt    = [ y[lag+1:end-1] ones(nobs-lag-1, 1) ]
        dfexo = 1
    elseif deterministic == :trend
        Δyt   = Δy[lag+1:end]
        xt    = [ y[lag+1:end-1] ones(nobs-lag-1, 1) (1:nobs-lag-1) ]
        dfexo = 2
    elseif deterministic == :squared_trend
        Δyt   = Δy[lag+1:end]
        xt    = [ y[lag+1:end-1] ones(nobs-lag-1, 1) (1:nobs-lag-1) [
                  i^2 for i = 1:nobs-lag-1] ]
        dfexo   = 3
    else
        throw(ArgumentError("deterministic = $(deterministic) is invalid"))
    end

    if lag > 0
        Δylag = hcat([Δy[lag-i+1:end-i] for i=1:lag]...)
        xt = [ xt Δylag ]
    end

    adf_coef = xt \ Δyt
    adf_res  = Δyt - xt * adf_coef
    n        = length(adf_res)
    k        = lag + 1 + dfexo
    sigma2   = dot(adf_res, adf_res) / (n - k)
    inv_xtxt = inv(xt' * xt)
    adf_stat = adf_coef[1] / sqrt(sigma2 * inv_xtxt[1,1])

    crit_val = adf_cv(nobs, deterministic)

    ADFTest(n, deterministic, lag, adf_stat, adf_coef[1], crit_val)

end

function adf_cv(nobs::Int, deterministic::Symbol)
    # computes critical values
    #
    # based on James G. MacKinnon, "Critical values for cointegration tests,"
    # QED Working Paper No. 1227, 2010. http://ideas.repec.org/p/qed/wpaper/1227.html

    cv_coeff = [
        # no constant, rows denote values for 1%, 5%, and 10% level, respectively
        -2.56574 -02.2358 -03.627    0.000
        -1.94100 -00.2686 -03.365   31.223
        -1.61682  00.2656 -02.714   25.364
        # constant
        -3.43035 -06.5393 -16.786 -079.433
        -2.86154 -02.8903 -04.234 -040.040
        -2.56677 -01.5384 -02.809    0.000
        # constant and trend
        -3.95877 -09.0531 -28.428 -134.155
        -3.41049 -04.3904 -09.036 -045.374
        -3.12705 -02.5856 -03.925 -022.380
        # constant, trend, and trend squared
        -4.37113 -11.5882 -35.819 -334.047
        -3.83239 -05.9057 -12.490 -118.284
        -3.55326 -03.6596 -05.293 -063.559
    ]

    # Equation A.1 in MacKinnon (2010)
    if deterministic == :none
        adf_crit_val = cv_coeff[1:3, 1] + cv_coeff[1:3, 2] ./ nobs +
                          cv_coeff[1:3, 3] ./ (nobs^2) + cv_coeff[1:3, 4] ./ (nobs^3)
    elseif deterministic == :constant
        adf_crit_val = cv_coeff[4:6, 1] + cv_coeff[4:6, 2] ./ nobs +
                          cv_coeff[4:6, 3] ./ (nobs^2) + cv_coeff[4:6, 4] ./ (nobs^3)
    elseif deterministic == :trend
        adf_crit_val = cv_coeff[7:9, 1] + cv_coeff[7:9, 2] ./ nobs +
                          cv_coeff[7:9, 3] ./ (nobs^2) + cv_coeff[7:9, 4] ./ (nobs^3)
    elseif deterministic == :squared_trend
        adf_crit_val = cv_coeff[10:12, 1] + cv_coeff[10:12, 2] ./ nobs +
                          cv_coeff[10:12, 3] ./ (nobs^2) + cv_coeff[10:12, 4] ./ (nobs^3)
    else
        throw(ArgumentError("deterministic = $(deterministic) is invalid"))
    end

end

function adf_pv_aux(adf_stat::Float64, deterministic::Symbol)
    # helper function for p-value computation
    #
    # based on James G. MacKinnon, "Approximate Asymptotic Distribution Functions for
    # Unit-Root and Cointegration Tests", Journal of Business & Economic Statistics,
    # Vol. 12, No. 2 (Apr., 1994), pp. 167-176, http://www.jstor.org/stable/1391481
    #
    # A slightly more accurate approach could be based on an algorithm described in
    # James G. MacKinnon, "Numerical Distribution Functions for Unit Root and Cointegration
    # Tests", Journal of Applied Econometrics, Vol. 11, No. 6 (Nov.-Dec., 1996),
    # pp. 601-618, http://www.jstor.org/stable/2285154
    # and made available under GPL as Fortran routines (used e.g. by R's urca package) at
    # http://qed.econ.queensu.ca/pub/faculty/mackinnon/numdist/

    pv_coeff_smallp = [
        0.6344 1.2378 0.032496 -1.04 -19.04
        2.1659 1.4412 0.038269 -1.61 -18.83
        3.2512 1.6047 0.049588 -2.89 -16.18
        4.0003 1.6580 0.048288 -3.21 -17.17
    ]

    pv_coeff_largep = [
        0.4797 0.93557 -0.06999  0.033066 -1.04 Inf
        1.7339 0.93202 -0.12745 -0.010368 -1.61 2.74
        2.5261 0.61654 -0.37956 -0.060285 -2.89 0.70
        3.0778 0.49529 -0.41477 -0.059359 -3.21 0.54
    ]

    if deterministic == :none
        tab_row = 1
    elseif deterministic == :constant
        tab_row = 2
    elseif deterministic == :trend
        tab_row = 3
    elseif deterministic == :squared_trend
        tab_row = 4
    else
        throw(ArgumentError("deterministic = $(deterministic) is invalid"))
    end

    if adf_stat < pv_coeff_smallp[tab_row, 5]
        aux_var = -Inf
    elseif adf_stat > pv_coeff_largep[tab_row, 6]
        aux_var = Inf
    else
        if adf_stat < pv_coeff_smallp[tab_row, 4]
            aux_var = pv_coeff_smallp[tab_row, 1] + pv_coeff_smallp[tab_row, 2] * adf_stat +
                        pv_coeff_smallp[tab_row, 3] * (adf_stat^2)
        else
            aux_var = pv_coeff_largep[tab_row, 1] + pv_coeff_largep[tab_row, 2] * adf_stat +
                        pv_coeff_largep[tab_row, 3] * (adf_stat^2) +
                        pv_coeff_largep[tab_row, 4] * (adf_stat^3)
        end

    end
end

testname(::ADFTest) = "Augmented Dickey-Fuller unit root test"
population_param_of_interest(x::ADFTest) =
    ("coefficient on lagged non-differenced variable", 0, x.coef)

function show_params(io::IO, x::ADFTest, ident)
    println(io, ident, "sample size in regression:          ", x.n)
    println(io, ident, "number of lags:                     ", x.lag)
    println(io, ident, "ADF statistic:                      ", x.stat)
    print(io, ident, "Critical values at 1%, 5%, and 10%: ")
    show(io, x.cv')
    println(io)
end

pvalue(x::ADFTest) = HypothesisTests.pvalue(Normal(0, 1), adf_pv_aux(x.stat, x.deterministic); tail=:left)
