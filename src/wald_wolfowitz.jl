export 
    WaldWolfowitzTest

struct WaldWolfowitzTest <: HypothesisTest
    N_above::Int
    N_below::Int
    N_runs::Int
    μ::Real
    σ::Real
    z::Real
end

testname(::WaldWolfowitzTest) = "Wald-Wolfowitz Test"
population_param_of_interest(x::WaldWolfowitzTest) = ("Number of runs", x.μ, x.N_runs) # parameter of interest: name, value under h0, point estimate
default_tail(::WaldWolfowitzTest) = :both


function show_params(io::IO, x::WaldWolfowitzTest, ident="")
    println(io, ident, "number of runs:  $(x.N_runs)")
    println(io, ident, "z-statistic:     $(x.z)")
end

"""
    WaldWolfowitzTest(x::AbstractVector{T<:Real})

Performs the Wald-Wolfowitz (or Runs) test to determine
randommness for a data-sequence. The null hypothesis is that the data is random, 
or independent.

Implements: [`pvalue`](@ref)
"""
function WaldWolfowitzTest(x::AbstractVector{T}) where T<:Real
    n = length(x)
    med = median(x)
    num_above = count(x .>= med)
    num_below = count(x .< med)

    # Get the expected value and standard deviation
    μ = 1 + 2 * num_above * num_below / n
    σ = sqrt((μ - 1) * (μ - 2) / (n - 1))

    # Get the number of runs
    signs = sign.(x .- med)
    num_runs = 1
    for k in 1:(n-1)
        if signs[k] != signs[k+1]
            num_runs += 1
        end
    end

    # calculate simple z-statistic
    z = (num_runs - μ) / σ
    WaldWolfowitzTest(num_above, num_below, num_runs, μ, σ, z)
end

pvalue(test::WaldWolfowitzTest; tail=:both) = pvalue(Normal(0.0, 1.0), test.z; tail=tail)

