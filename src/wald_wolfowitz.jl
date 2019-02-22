export 
    WaldWolfowitzTest
    TwoValuedWaldWolfowitzTest

struct WaldWolfowitzTest <: HypothesisTest
    nabove::Int
    nbelow::Int
    nruns::Int
    μ::Real
    σ::Real
    z::Real
end

testname(::WaldWolfowitzTest) = "Wald-Wolfowitz Test"
population_param_of_interest(x::WaldWolfowitzTest) = ("Number of runs", x.μ, x.nruns) # parameter of interest: name, value under h0, point estimate
default_tail(::WaldWolfowitzTest) = :both


function show_params(io::IO, x::WaldWolfowitzTest, ident="")
    println(io, ident, "number of runs:  $(x.nruns)")
    println(io, ident, "z-statistic:     $(x.z)")
end

"""
    WaldWolfowitzTest(x::AbstractVector{<:Real})

Performs the Wald-Wolfowitz (or Runs) test of the null hypothesis that the given data is random, or independent.
The data is transformed to two-valued data by labelling a point as above or below the median of the data.

Implements: [`pvalue`](@ref)
"""
function WaldWolfowitzTest(x::AbstractVector{T}) where T<:Real
    med = median(x)
    transformed = x .>= med
    TwoValuedWaldWolfowitzTest(transformed)
end

"""
    TwoValuedWaldWolfowitzTest(x::Vector{::bool})

Performs the Wald-Wolfowitz (or Runs) test of the null hypothesis that the given data is random, or independent.
This test assumes two-valued boolean data.

Implements: [`pvalue`](@ref)
"""
function TwoValuedWaldWolfowitzTest(x::Vector{::bool))
    n = length(x)
    nabove = sum(x)
    nbelow = n - nabove
    
    # Get the expected value and standard deviation
    μ = 1 + 2 * nabove * (nbelow / n)
    σ = sqrt((μ - 1) * (μ - 2) / (n - 1))

    # Get the number of runs
    nruns = 1
    for k in 1:(n-1)
        if x[k] != x[k+1]
            nruns += 1
        end
    end

    # calculate simple z-statistic
    z = (nruns - μ) / σ
    WaldWolfowitzTest(nabove, nbelow, nruns, μ, σ, z)

end

pvalue(test::WaldWolfowitzTest; tail=:both) = pvalue(Normal(0.0, 1.0), test.z; tail=tail)

