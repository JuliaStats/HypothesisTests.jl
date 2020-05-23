# Tests for Equality of Variances

export LeveneTest, BrownForsytheTest, FlignerKilleenTest

abstract type VarianceEqualityTest <: HypothesisTest end

population_param_of_interest(::VarianceEqualityTest) =
    ("Equality of variances", NaN, NaN)

# Levene Test
struct LeveneTest <: VarianceEqualityTest
    Nᵢ::Vector{Int}         # number of observations in each group
    W::Float64              # test statistic: F statistic
end

"""
    LeveneTest(groups::AbstractVector{<:Real}...; statistic=mean)

Perform Levene's test of the hypothesis that that the `groups` variances are equal.
"""
function LeveneTest(groups::AbstractVector{<:Real}...; statistic=mean)
    Nᵢ = [length(g) for g in groups]
    N = sum(Nᵢ)
    k = length(Nᵢ)
    Ȳ = [statistic(g) for g in groups]
    Zᵢⱼ = [abs.(g .- Ȳ[i]) for (i,g) in enumerate(groups)]
    Zᵢ₋ = map(mean, Zᵢⱼ)
    Z₋₋ = mean(Zᵢ₋)
    SSgrp = sum(Nᵢ .* (Zᵢ₋ .- Z₋₋).^2)
    SSerr = sum( sum.((Z .- Zᵢ₋[i]).^2 for (i,Z) in enumerate(Zᵢⱼ)) )
    W = ((N-k)*SSgrp)/((k-1)*SSerr)
    LeveneTest(Nᵢ,W)
end

"""
    BrownForsytheTest(groups::AbstractVector{<:Real}...)

The Brown–Forsythe test is a statistical test for the equality of `groups` variances.

*Note: This is essentially, the Levene's test with the median instead of the mean statistic for computing the spread within each group.*
"""
BrownForsytheTest(groups::AbstractVector{<:Real}...) = LeveneTest(groups...; statistic=median)

StatsBase.nobs(L::LeveneTest) = L.Nᵢ
StatsBase.dof(L::LeveneTest) = let k = length(L.Nᵢ); (k-1, sum(nobs(L))-k) end

testname(::LeveneTest) = "Levene's Test for Equality of Variances"
default_tail(::LeveneTest) = :right
pvalue(L::LeveneTest; tail=:right) = pvalue(FDist(dof(L)...), L.W, tail=tail)

function show_params(io::IO, L::LeveneTest, indent="")
    println(io, indent, "number of observations: ", nobs(L))
    println(io, indent, "W statistic:            ", L.W)
    println(io, indent, "degrees of freedom:     ", dof(L))
end

# Fligner-Killeen Test
struct FlignerKilleenTest <: VarianceEqualityTest
    Nᵢ::Vector{Int}         # number of observations in each group
    FK::Float64             # test statistic: chi-square statistic
end

"""
    FlignerKilleenTest(groups::AbstractVector{<:Real}...)

Perform Fligner-Killeen median test of the null hypothesis that the `groups`
have equal variances, a test for homogeneity of variances.

It is a %k%-sample simple linear rank method that uses the ranks of the absolute values of the centered samples, and weights
```math
a_{N,i} = \\Phi^{-1}(1/2 + (i/2(N+1)))
```

Implements: [`pvalue`](@ref)

# References

  * Conover, W. J., Johnson, M. E., Johnson, M. M., A comparative study of tests
    for homogeneity of variances, with applications to the outer continental shelf bidding data.
    Technometrics, 23, 351–361, 1980

# External links

  * [Fligner-Killeen test on Statistical Analysis Handbook
    ](https://www.statsref.com/HTML/index.html?fligner-killeen_test.html)
"""
function FlignerKilleenTest(groups::AbstractVector{<:Real}...)
    Nᵢ = [length(g) for g in groups]
    N = sum(Nᵢ)
    k = length(Nᵢ)
    Ȳ = [median(g) for g in groups]
    Zᵢⱼ = [abs.(g .- Ȳ[i]) for (i,g) in enumerate(groups)]
    (ranks, tieadj) = tiedrank_adj(vcat(Zᵢⱼ...))
    qᵢⱼ = quantile.(Normal(),0.5 .+ ranks./2(N+1))
    qᵢ₋ = vec(mean(reshape(qᵢⱼ, :, k), dims=1))
    WSS = Nᵢ .* (qᵢ₋ .- mean(qᵢ₋)).^2
    FK = sum(WSS) / var(qᵢⱼ)
    FlignerKilleenTest(Nᵢ, FK)
end

StatsBase.nobs(t::FlignerKilleenTest) = t.Nᵢ
StatsBase.dof(t::FlignerKilleenTest) = length(t.Nᵢ)-1

testname(::FlignerKilleenTest) = "Fligner-Killeen median test for homogeneity of variances"
default_tail(::FlignerKilleenTest) = :right
pvalue(t::FlignerKilleenTest; tail=:right) = pvalue(Chisq(dof(t)), t.FK, tail=tail)

function show_params(io::IO, t::FlignerKilleenTest, indent="")
    println(io, indent, "number of observations: ", nobs(t))
    println(io, indent, "FK statistic:           ", t.FK)
    println(io, indent, "degrees of freedom:     ", dof(t))
end
