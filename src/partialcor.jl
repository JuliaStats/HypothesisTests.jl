# Partial correlation

export PartialCorTest

"""
    PartialCorTest(x, y, Z)

Perform a t-test for the hypothesis that ``\\text{Cor}(x,y|Z=z) = 0``, i.e. the partial
correlation of vectors `x` and `y` given the matrix `Z` is zero.

Implements `pvalue` and `confint`. See also [`partialcor`](@ref).
"""
struct PartialCorTest <: HypothesisTest
    r::Real
    n::Int
    k::Int
    t::Real

    # Error checking is done in `partialcor`
    function PartialCorTest(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix)
        r = partialcor(x, y, Z)
        n, k = size(Z)
        t = r * sqrt((n - 2 - k) / (1 - r^2))
        return new(r, n, k, t)
    end

    function PartialCorTest(x::AbstractVector, y::AbstractVector, z::AbstractVector)
        r = partialcor(x, y, z)
        n = length(z)
        t = r * sqrt((n - 3) / (1 - r^2))
        return new(r, n, 1, t)
    end
end

testname(::PartialCorTest) = "Test for partial correlation"
population_param_of_interest(p::PartialCorTest) = ("Partial correlation", zero(p.r), p.r)

StatsBase.nobs(p::PartialCorTest) = p.n
StatsBase.dof(p::PartialCorTest) = p.n - 2 - p.k

function StatsBase.confint(test::PartialCorTest, alpha::Float64=0.05)
    dof(test) > 1 || return (-one(test.r), one(test.r)) # Otherwise we can get NaNs
    q = quantile(Normal(), 1 - alpha / 2)
    fisher = (log1p(test.r) - log1p(-test.r)) / 2
    bound = q / sqrt(test.n - 3 - test.k)
    lo = 2 * (fisher - bound)
    hi = 2 * (fisher + bound)
    elo = Statistics.clampcor(expm1(lo) / (exp(lo) + 1))
    ehi = Statistics.clampcor(expm1(hi) / (exp(hi) + 1))
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
