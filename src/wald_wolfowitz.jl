export WaldWolfowitzTest

struct WaldWolfowitzTest{T <: Real} <: HypothesisTest
    nabove::Int     # Number of points above median (or of value a)
    nbelow::Int     # Number of points below median (or of value b)
    nruns::Int      # Number of runs
    μ::T            # Expected mean
    σ::T            # Expected variance
    z::T            # test z-statistic from Normal distribution
end

testname(::WaldWolfowitzTest) = "Wald-Wolfowitz Test"
population_param_of_interest(x::WaldWolfowitzTest) = ("Number of runs", x.μ, x.nruns) # parameter of interest: name, value under h0, point estimate
default_tail(::WaldWolfowitzTest) = :both
pvalue(test::WaldWolfowitzTest; tail=:both) = pvalue(Normal(0.0, 1.0), test.z; tail=tail)


function show_params(io::IO, x::WaldWolfowitzTest, ident="")
    println(io, ident, "number of runs:  $(x.nruns)")
    println(io, ident, "z-statistic:     $(x.z)")
end

"""
    WaldWolfowitzTest(x::AbstractVector{Bool})
    WaldWolfowitzTest(x::AbstractVector{<:Real})

Perform the Wald-Wolfowitz (or Runs) test of the null hypothesis that the given data is random, or independently sampled.
The data can come as many-valued or two-valued (Boolean). If many-valued, the sample is transformed by labelling each
element as above or below the median.

Implements: [`pvalue`](@ref)
"""
function WaldWolfowitzTest(x::AbstractVector{Bool})
    n = length(x)
    nabove = sum(x)
    nbelow = n - nabove
    
    # Get the expected value and standard deviation
    μ = 1 + 2 * nabove * (nbelow / n)
    σ = sqrt((μ - 1) * (μ - 2) / (n - 1))

    # Get the number of runs
    nruns = 1
    for k in 1:(n - 1)
        @inbounds if x[k] != x[k + 1]
            nruns += 1
        end
    end

    # calculate simple z-statistic
    z = (nruns - μ) / σ
    WaldWolfowitzTest(nabove, nbelow, nruns, μ, σ, z)
end

WaldWolfowitzTest(x::AbstractVector{<:Real}) = WaldWolfowitzTest(x .>= median(x))
