# Partial correlation

export PartialCorTest

"""
    PartialCorTest(x, y, Z)

Perform a t-test for the hypothesis that ``\\text{Cor}(x,y|Z=z) = 0``, i.e. the partial
correlation of vectors `x` and `y` given the matrix `Z` is zero.

Implements `pvalue` and `confint`. See also [`partialcor`](@ref).
"""
struct PartialCorTest{T<:Real} <: HypothesisTest
    r::T
    n::Int
    k::Int
    t::T

    # Error checking is done in `partialcor`
    function PartialCorTest(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix)
        r = partialcor(x, y, Z)
        n, k = size(Z)
        t = r * sqrt((n - 2 - k) / (1 - r^2))
        return new{typeof(r)}(r, n, k, t)
    end

    function PartialCorTest(x::AbstractVector, y::AbstractVector, z::AbstractVector)
        r = partialcor(x, y, z)
        n = length(z)
        t = r * sqrt((n - 3) / (1 - r^2))
        return new{typeof(r)}(r, n, 1, t)
    end
end

testname(::PartialCorTest) = "Test for partial correlation"
population_param_of_interest(p::PartialCorTest) = ("Partial correlation", zero(p.r), p.r)

StatsBase.nobs(p::PartialCorTest) = p.n
StatsBase.dof(p::PartialCorTest) = p.n - 2 - p.k

function StatsBase.confint(test::PartialCorTest{T}, alpha::Float64=0.05) where T
    df = dof(test)
    df > 1 || return (-one(T), one(T))  # Otherwise we can get NaNs
    q = quantile(TDist(df), 1 - alpha / 2)
    fisher = atanh(test.r)
    bound = q / sqrt(test.n - 3 - test.k)
    lo = 2 * (fisher - bound)
    hi = 2 * (fisher + bound)
    elo = Statistics.clampcor(tanh(fisher - bound))
    ehi = Statistics.clampcor(tanh(fisher + bound))
    return (elo, ehi)
end

default_tail(::PartialCorTest) = :both
pvalue(test::PartialCorTest; tail=:both) = pvalue(TDist(dof(test)), test.t, tail=tail)

function show_params(io::IO, test::PartialCorTest, indent="")
    println(io, indent, "number of observations:          ", nobs(test))
    println(io, indent, "number of conditional variables: ", test.k)
    println(io, indent, "t-statistic:                     ", test.t)
    println(io, indent, "degrees of freedom:              ", dof(test))
end
