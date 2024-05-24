# KPSS.jl
# Kwiatkowski–Phillips–Schmidt–Shin unit root test
#
# Copyright (C) 2024 Jared Schwartz
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

export KPSSTest

struct KPSSTest <: HypothesisTest
    stat::Float64  # test statistic
    regression::Symbol  # regression type
    lag::Int   # number of lags  
end

"""
    KPSSTest(y::AbstractVector{d}, regression::Symbol, lag::Union{Symbol, Int}) where d<:Real

Compute the Kwiatkowski-Phillips-Schmidt-Shin unit root test.

`y` is the time series to be tested, `regression` indicates constant or constant+trend
(options: `:constant`, `:trend`) and `lag` is the number of lagged
first-differences included in the test regression. `lag` can also be configured to automatically find 
the find the optimal number of lags for the data using `:auto` for the the Hobijn et al. data-dependent method
or `:legacy` for the the KPSS data-independent method.

# References

  * Hobijn, Bart, Philip Hans Franses, and Marius Ooms. "Generalizations of the KPSS-Test for Stationarity." 
    Statistica Neerlandica 58, no. 4 (2004): 483-502. https://doi.org/10.1111/j.1467-9574.2004.00272.x.

  * Kwiatkowski, Denis, Peter C. B. Phillips, Peter Schmidt, and Yongcheol Shin. "Testing the Null Hypothesis of
    Stationarity against the Alternative of a Unit Root: How Sure Are We That Economic Time Series Have a Unit Root?" 
    Journal of Econometrics 54, no. 1 (October 1, 1992): 159-78. https://doi.org/10.1016/0304-4076(92)90104-Y.
  
# External links

  * [KPSS test on Wikipedia](https://en.wikipedia.org/wiki/KPSS_test)
"""

function autolag(resids::AbstractVector{d}, T::Int) where d <: Real
    # Computes the optimal number of lags using the Bartlett kernel as described in Hobijn et al. (2004).
    # Reference: Table 3 p.489 of Hobijn et al. (2004)

    n = trunc(Int, T^(2.0 / 9.0))
    
    ŝ_0 = sum(resids .^ 2) / T
    ŝ_j = 0.0
    for i in 1:n
        resids_prod = dot(resids[i+1:end], resids[1:T-i])
        resids_prod /= T / 2.0
        ŝ_0 += resids_prod
        ŝ_j += i * resids_prod
    end

    γ̂ = 1.1447 * ((ŝ_j / ŝ_0)^2)^(1.0 / 3.0)
    m̂ₜ = min(T, trunc(Int,γ̂ * T^(1.0 / 3.0)))
    return m̂ₜ

end

function s²_l(ϵ::AbstractVector{d}, T::Int, l::Int) where d <: Real
    # calculate weighted s² using Bartlett kernel weighting as in Kwiatkowski et al. (1992) and Hobijn et al. (2004)
    # Reference: Equation (10), p. 164 (Kwiatkowski et al., 1992).

    SSE = sum(ϵ.^2)

    weighted_sum = 0
    for s in 1:l
        weight = 1.0 - (s / (l + 1.0)) 
        resids_prod = dot(ϵ[s+1:end], ϵ[1:T-s])
        weighted_sum += 2 * weight * resids_prod
    end
    
    ŝ² = (SSE + weighted_sum) / T
    
    return ŝ²
end


function KPSSTest(y::AbstractVector{d}, regression::Symbol= :constant, lag::Union{Symbol, Int}= :auto) where d <: Real
    
    # Number of observations
    T = length(y)
    
    # Handle regression argument and detrend based on regression type
    if regression == :constant
        ϵ = y .- mean(y)
    elseif regression == :trend
        t = collect(1:T)
        X = [ones(T) t]
        β = (X'X)\(X'y)
        ϵ = y .- X*β
    else
        throw(ArgumentError("regression = $(regression) is invalid. Choose :constant for constant or :trend for constant and trend."))
    end

    # Handle lags argument
    if lag == :auto
        # data dependent method defined in Hobijn et al. (2004)
        lag = autolag(y, T)
    elseif lag == :legacy
        # data independent method used for Kwiatkowski et al. (1992)
        lag = trunc(Int, ceil(12.0 * (T / 100.0)^(1 / 4.0)))  
    elseif typeof(lag) == Int
        lag = lag
    else
        throw(ArgumentError("lag = $(lag) is invalid. lag must be :legacy, :auto, or an Int"))
    end

    lag = min(lag, T) # lags cannot be greater than the length of the series

    η  = sum(cumsum(ϵ).^2) / (T^2)
    ŝ² = s²_l(ϵ, T, lag)

    stat = η / ŝ²
    return KPSSTest(stat, regression, lag)
end


function interp_p(kpss_stat::Real, regression::Symbol)
    
    # Approximate p-values in Table 1 p.166 (Kwiatkowski et al., 1992)
    p_values = [0.10, 0.05, 0.025, 0.01]
    critical_values = Dict(
        :constant => [0.347, 0.463, 0.574, 0.739],
        :trend => [0.119, 0.146, 0.176, 0.216]
    )

    if kpss_stat > maximum(critical_values[regression])
        @warn "p-value is outside of interpolation table. It is smaller than the returned p-value"
        return 0.01
    elseif kpss_stat < minimum(critical_values[regression])
        @warn "p-value is outside of interpolation table. It is greater than the returned p-value"
        return 0.10
    else
       # Manual linear interpolation
       x = kpss_stat
       xp = critical_values[regression]
       fp = p_values
       function interp(x::d, xp::Vector{d}, fp::Vector{d}) where d <: Real
            for i in 1:length(xp)-1
                if xp[i] <= x <= xp[i+1]
                    p_value =  fp[i] + (fp[i+1] - fp[i]) * (x - xp[i]) / (xp[i+1] - xp[i])
                end
            end
        
            return p_value
        end

        p_value = interp(x, xp, fp)
    end
    return p_value
end

StatsAPI.pvalue(x::KPSSTest) = interp_p(x.stat, x.regression)

testname(::KPSSTest) = "Kwiatkowski-Phillips-Schmidt-Shin unit root test"

function show_params(io::IO, x::KPSSTest, ident)
    println(io, ident, "Number of lags:                     ", x.lag)
    println(io, ident, "Regression type:                    ", x.regression)
    println(io, ident, "KPSS statistic:                     ", x.stat)
    println(io)
end