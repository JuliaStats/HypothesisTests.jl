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
By default the mean `statistic` is used for centering in each of the `groups`, but
other statistics are accepted: median or truncated mean, see [`BrownForsytheTest`](@ref).

The test statistic, ``W``, is equivalent to the ``F`` statistic, and is defined as follows:

```math
W = \\frac{(N-k)}{(k-1)} \\cdot \\frac{\\sum_{i=1}^k N_i (Z_{i\\cdot}-Z_{\\cdot\\cdot})^2} {\\sum_{i=1}^k \\sum_{j=1}^{N_i} (Z_{ij}-Z_{i\\cdot})^2},
```
where

* ``k`` is the number of different groups to which the sampled cases belong,
* ``N_i`` is the number of cases in the ``i``th group,
* ``N`` is the total number of cases in all groups,
* ``Y_{ij}`` is the value of the measured variable for the ``j``th case from the ``i``th group,
* ``Z_{ij} = |Y_{ij} - \\bar{Y}_{i\\cdot}|``, ``\\bar{Y}_{i\\cdot}`` is a mean of the  ``i``th group,
* ``Z_{i\\cdot} = \\frac{1}{N_i} \\sum_{j=1}^{N_i} Z_{ij}`` is the mean of the ``Z_{ij}`` for group ``i``,
* ``Z_{\\cdot\\cdot} = \\frac{1}{N} \\sum_{i=1}^k \\sum_{j=1}^{N_i} Z_{ij}`` is the mean of all ``Z_{ij}``.

The test statistic ``W`` is approximately ``F``-distributed with ``k-1`` and ``N-k`` degrees of freedom.

# References

  * Levene, Howard, "Robust tests for equality of variances". In Ingram Olkin; Harold Hotelling; et al. (eds.).
     Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling.
     Stanford University Press. pp. 278–292, 1960

# External links

  * [Levene's test on Wikipedia
    ](https://en.wikipedia.org/wiki/Levene%27s_test)
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

The Brown–Forsythe test is a modification of the Levene's test with the median instead of the mean statistic for computing the spread within each group.

Implements: [`pvalue`](@ref)

# References

  * Brown, Morton B.; Forsythe, Alan B., "Robust tests for the equality of variances".
    Journal of the American Statistical Association. 69: 364–367, 1974
    doi:[10.1080/01621459.1974.10482955](https://doi.org/10.1080%2F01621459.1974.10482955).

# External links

  * [Brown–Forsythe test on Wikipedia
    ](https://en.wikipedia.org/wiki/Brown%E2%80%93Forsythe_test)
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

This test is most robust against departures from normality, see references.
It is a ``k``-sample simple linear rank method that uses the ranks of the absolute values of the centered samples and weights
```math
a_{N,i} = \\Phi^{-1}(1/2 + (i/2(N+1)))
```
The version implemented here uses median centering in each of the samples.

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
