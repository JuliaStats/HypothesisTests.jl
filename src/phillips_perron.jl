# phillips_perron.jl
# Phillips-Perron unit root test

export PhillipsPerronTest

struct PhillipsPerronTest <: HypothesisTest
    n               ::Int               # number of observations in the test regression
    deterministic   ::Symbol            # deterministic terms included in regression
    lag             ::Int               # number of lags in the long-run variance estimate
    stat            ::Float64           # test statistic Z_t
    coef            ::Float64           # coefficient on lagged (non-differenced) variable
    cv              ::Vector{Float64}   # critical values
end

# The "short" lag rule of Schwert (1989), the default of R's urca::ur.pp and tseries::pp.test
pp_default_lag(n::Integer) = n > 0 ? floor(Int, 4 * (n / 100)^(1 / 4)) : 0

"""
    PhillipsPerronTest(y::AbstractVector{<:Real}; deterministic::Symbol=:constant,
                       lag::Integer=floor(Int, 4 * ((length(y) - 1) / 100)^(1 / 4)))

Compute the Phillips-Perron unit root test.

`y` is the time series to be tested. The null hypothesis is that `y` has a unit root. The
test estimates the Dickey-Fuller regression without lagged differences,
```math
\\Delta y_t = (\\rho - 1) y_{t-1} + d_t + u_t, \\qquad t = 2, \\ldots, T,
```
where the deterministic terms ``d_t`` are set by `deterministic` (options: `:none`,
`:constant`, `:trend`). Where [`ADFTest`](@ref) adds lagged differences to the regression
to absorb serial correlation in ``u_t``, this test keeps the regression as it is and
corrects the t-statistic of ``\\hat\\rho - 1`` instead:
```math
Z_t = \\sqrt{\\frac{\\hat\\gamma_0}{\\hat\\lambda^2}} \\frac{\\hat\\rho - 1}{\\hat\\sigma_{\\hat\\rho}}
    - \\frac{\\hat\\lambda^2 - \\hat\\gamma_0}{2 \\hat\\lambda} \\frac{n \\hat\\sigma_{\\hat\\rho}}{s},
```
where ``n = T - 1`` is the number of observations in the regression,
``\\hat\\sigma_{\\hat\\rho}`` is the OLS standard error of ``\\hat\\rho``, ``s^2`` is the OLS
residual variance (with ``n - k`` degrees of freedom for ``k`` regressors),
``\\hat\\gamma_0 = n^{-1} \\sum_t \\hat u_t^2``, and ``\\hat\\lambda^2`` is the Newey-West
estimate of the long-run variance of ``u_t`` with Bartlett weights and `lag` lags,
```math
\\hat\\lambda^2 = \\hat\\gamma_0
    + 2 \\sum_{j=1}^{\\mathrm{lag}} \\left(1 - \\frac{j}{\\mathrm{lag} + 1}\\right) \\hat\\gamma_j,
\\qquad
\\hat\\gamma_j = \\frac{1}{n} \\sum_{t=j+1}^{n} \\hat u_t \\hat u_{t-j}.
```
This is the form of the statistic given by Hamilton (1994). With `lag = 0`, ``Z_t`` is the
Dickey-Fuller t-statistic.

The default `lag` is ``\\lfloor 4 (n/100)^{1/4} \\rfloor``, the short rule of Schwert (1989)
that R's `urca::ur.pp` and `tseries::pp.test` use by default. Python's
`arch.unitroot.PhillipsPerron` uses ``\\lceil 12 (T/100)^{1/4} \\rceil`` lags by default.

``Z_t`` has the same asymptotic distribution as the Dickey-Fuller t-statistic, so the
critical values and p-values are the ones of [`ADFTest`](@ref): response surface
regressions following MacKinnon (2010), evaluated at ``n``, and MacKinnon (1994),
respectively.

# References

  * [Phillips and Perron (1988). Testing for a unit root in time series
    regression](@cite phillips1988)
  * [Newey and West (1987). A Simple, Positive Semi-Definite, Heteroskedasticity and
    Autocorrelation Consistent Covariance Matrix](@cite newey1987)
  * [Hamilton (1994). Time Series Analysis](@cite hamilton1994)
  * [Schwert (1989). Tests for Unit Roots: A Monte Carlo Investigation](@cite schwert1989)
  * [MacKinnon (2010). Critical values for cointegration tests](@cite mackinnon2010)
  * [MacKinnon (1994). Approximate Asymptotic Distribution Functions for
    Unit-Root and Cointegration Tests](@cite mackinnon1994)

# External links

  * [Phillips-Perron test on Wikipedia](https://en.wikipedia.org/wiki/Phillips–Perron_test)
"""
function PhillipsPerronTest(y::AbstractVector{<:Real}; deterministic::Symbol=:constant,
                            lag::Integer=pp_default_lag(length(y) - 1))
    n = length(y) - 1
    Δy = Vector{Float64}(diff(y))
    ylag = Vector{Float64}(y[1:end-1])

    if deterministic === :none
        X = reshape(ylag, n, 1)
    elseif deterministic === :constant
        X = [ylag ones(n)]
    elseif deterministic === :trend
        X = [ylag ones(n) 1:n]
    else
        throw(ArgumentError("deterministic = $(deterministic) is invalid"))
    end
    k = size(X, 2)
    n > k || throw(ArgumentError("y must have at least $(k + 2) observations"))
    0 <= lag < n || throw(ArgumentError("lag = $(lag) is invalid, it must be between 0 and $(n - 1)"))

    F = qr(X)
    coef = F \ Δy
    res = Δy - X * coef
    γ0 = dot(res, res) / n
    s2 = dot(res, res) / (n - k)

    # Newey-West estimate of the long-run variance of the errors (Bartlett kernel)
    λ2 = γ0
    for j in 1:lag
        λ2 += 2 * (1 - j / (lag + 1)) * dot(view(res, j+1:n), view(res, 1:n-j)) / n
    end

    # first diagonal element of inv(X'X) = inv(R) * inv(R)'
    v11 = sum(abs2, view(inv(UpperTriangular(F.R)), 1, :))
    se = sqrt(s2 * v11)
    # n * se / s simplifies to n * sqrt(v11)
    stat = sqrt(γ0 / λ2) * coef[1] / se - (λ2 - γ0) / (2 * sqrt(λ2)) * n * sqrt(v11)

    PhillipsPerronTest(n, deterministic, lag, stat, coef[1], adf_cv(n, deterministic))
end

testname(::PhillipsPerronTest) = "Phillips-Perron unit root test"
population_param_of_interest(x::PhillipsPerronTest) =
    ("coefficient on lagged non-differenced variable", 0, x.coef)
default_tail(::PhillipsPerronTest) = :left

function show_params(io::IO, x::PhillipsPerronTest, ident)
    println(io, ident, "sample size in regression:          ", x.n)
    println(io, ident, "lags in long-run variance:          ", x.lag)
    println(io, ident, "Z_t statistic:                      ", x.stat)
    print(io, ident, "Critical values at 1%, 5%, and 10%: ")
    join(io, x.cv, ", ")
    println(io)
end

StatsAPI.pvalue(x::PhillipsPerronTest) =
    pvalue(Normal(0, 1), adf_pv_aux(x.stat, x.deterministic); tail=:left)
