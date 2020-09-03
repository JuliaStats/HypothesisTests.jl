# Tests for Equality of Variances

export OneWayANOVATest, LeveneTest, BrownForsytheTest, FlignerKilleenTest

struct VarianceEqualityTest{TD <: ContinuousUnivariateDistribution} <: HypothesisTest
    Nᵢ::Vector{Int}
    SStᵢ::Vector{Float64}
    SSeᵢ::Vector{Float64}
    DFt::Int
    DFe::Int
    # test name, parameter of interest, test statistic name
    description::Tuple{String, String, String}
end

population_param_of_interest(t::VarianceEqualityTest) =
    (t.description[2], "all equal", NaN)
testname(t::VarianceEqualityTest) = t.description[1]
teststatisticname(t::VarianceEqualityTest{TD}) where {TD <: ContinuousDistribution} =
    length(t.description[3]) != 0 ? t.description[3] : (TD <: FDist ? "F" : "χ²")

StatsBase.nobs(t::VarianceEqualityTest) = t.Nᵢ
StatsBase.dof(t::VarianceEqualityTest{Chisq}) = t.DFt
StatsBase.dof(t::VarianceEqualityTest{FDist}) = (t.DFt, t.DFe)

teststatistic(t::VarianceEqualityTest{FDist}) = (t.DFe/t.DFt)*sum(t.SStᵢ)/sum(t.SSeᵢ)
function teststatistic(t::VarianceEqualityTest{Chisq})
    y = sum(t.SStᵢ)/sum(t.SSeᵢ)
    y*(t.DFe+t.DFt)/(1 + y) # sum(t.SStᵢ)/t.s²
end
pvalue(t::VarianceEqualityTest{TD}; tail=:right) where {TD <: ContinuousDistribution} =
    pvalue(TD(dof(t)...), teststatistic(t), tail=tail)

function show_params(io::IO, t::VarianceEqualityTest{TD},
                     indent="") where {TD <: ContinuousDistribution}
    println(io, indent, "number of observations: ", nobs(t))
    println(io, indent, rpad("$(teststatisticname(t)) statistic:", 24), teststatistic(t))
    println(io, indent, "degrees of freedom:     ", dof(t))
end

function Base.show(io::IOContext, t::VarianceEqualityTest)
    if !get(io, :table, false) # No table
        show(io.io, t)
    else
        println(io, testname(t))
        println(io, repeat("-", 55))
        SSt = sum(t.SStᵢ)
        SSe = sum(t.SSeᵢ)
        MSt = SSt/t.DFt
        MSe = SSe/t.DFe
        println(io, "Source            SS    DF        MS         F  P-value")
        println(io, repeat("-", 55))
        StatsBase.@printf(io, "Treatments  %8.3f  %4d  %8.3f  %8.3f  %7.5f\n",
                              SSt, t.DFt, MSt, MSt/MSe, pvalue(t))
        StatsBase.@printf(io, "Error       %8.3f  %4d  %8.3f\n", SSe, t.DFe, MSe)
        println(io, repeat("-", 55))
        StatsBase.@printf(io, "Total       %8.3f  %4d\n", SSt+SSe, t.DFt+t.DFe)
    end
end

function anova(scores::AbstractVector{<:Real}...)
    Nᵢ = [length(g) for g in scores]
    Z̄ᵢ = mean.(scores)
    Z̄ = mean(Z̄ᵢ)
    SStᵢ = Nᵢ .* (Z̄ᵢ .- Z̄).^2
    SSeᵢ = sum.( (z .- z̄).^2 for (z, z̄) in zip(scores, Z̄ᵢ) )
    (Nᵢ, SStᵢ, SSeᵢ)
end

"""
    OneWayANOVATest(groups::AbstractVector{<:Real}...)

Perform one-way analysis of variance test of the hypothesis that that the `groups`
means are equal.

The one-way analysis of variance (one-way ANOVA) is a technique that can be used
to compare means of two or more samples. The ANOVA tests the null hypothesis, which states
that samples in all groups are drawn from populations with the same mean values.
To do this, two estimates are made of the population variance.
The ANOVA produces an F-statistic, the ratio of the variance calculated among
the means to the variance within the samples.

Implements: [`pvalue`](@ref)

# External links

  * [One-way analysis of variance on Wikipedia
    ](https://en.wikipedia.org/wiki/One-way_analysis_of_variance)
"""
function OneWayANOVATest(groups::AbstractVector{<:Real}...)
    Nᵢ, SStᵢ, SSeᵢ = anova(groups...)
    k = length(Nᵢ)
    VarianceEqualityTest{FDist}(Nᵢ, SStᵢ, SSeᵢ, k-1, sum(Nᵢ)-k,
        ("One-way analysis of variance (ANOVA) test","Means","F"))
end

"""
    LeveneTest(groups::AbstractVector{<:Real}...; scorediff=abs, statistic=mean)

Perform Levene's test of the hypothesis that that the `groups` variances are equal.
By default the mean `statistic` is used for centering in each of the `groups`, but
other statistics are accepted: median or truncated mean, see [`BrownForsytheTest`](@ref).
By default the absolute value of the score difference, `scorediff`, is used, but
other functions are accepted: x² or √|x|.

The test statistic, ``W``, is equivalent to the ``F`` statistic, and is defined as follows:

```math
W = \\frac{(N-k)}{(k-1)} \\cdot \\frac{\\sum_{i=1}^k N_i (Z_{i\\cdot}-Z_{\\cdot\\cdot})^2}
    {\\sum_{i=1}^k \\sum_{j=1}^{N_i} (Z_{ij}-Z_{i\\cdot})^2},
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

Implements: [`pvalue`](@ref)

# References

  * Levene, Howard, "Robust tests for equality of variances".
     In Ingram Olkin; Harold Hotelling; et al. (eds.).
     Contributions to Probability and Statistics: Essays in Honor of Harold Hotelling.
     Stanford University Press. pp. 278–292, 1960

# External links

  * [Levene's test on Wikipedia
    ](https://en.wikipedia.org/wiki/Levene%27s_test)
"""
function LeveneTest(groups::AbstractVector{<:Real}...; scorediff=abs, statistic=mean)
    # calculate scores
    Zᵢⱼ = [scorediff.(g .- statistic(g)) for g in groups]
    # anova
    Nᵢ, SStᵢ, SSeᵢ = anova(Zᵢⱼ...)
    k = length(Nᵢ)
    VarianceEqualityTest{FDist}(Nᵢ, SStᵢ, SSeᵢ, k-1, sum(Nᵢ)-k,
        ("Levene's test","Variances","W"))
end

"""
    BrownForsytheTest(groups::AbstractVector{<:Real}...)

The Brown–Forsythe test is a statistical test for the equality of `groups` variances.

The Brown–Forsythe test is a modification of the Levene's test with the median instead
of the mean statistic for computing the spread within each group.

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

"""
    FlignerKilleenTest(groups::AbstractVector{<:Real}...)

Perform Fligner-Killeen median test of the null hypothesis that the `groups`
have equal variances, a test for homogeneity of variances.

This test is most robust against departures from normality, see references.
It is a ``k``-sample simple linear rank method that uses the ranks of the absolute values
of the centered samples and weights
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
    # calculate scores
    Zᵢⱼ = [abs.(g .- median(g)) for g in groups]
    # rank scores
    (ranks, tieadj) = tiedrank_adj(vcat(Zᵢⱼ...))
    qᵢⱼ = quantile.(Normal(),0.5 .+ ranks./2(length(ranks)+1))
    Nᵢ = pushfirst!(cumsum([length(g) for g in groups]),0)
    Qᵢⱼ = [qᵢⱼ[(Nᵢ[i]+1):(Nᵢ[i+1])] for i in 1:length(Nᵢ)-1]
    # anova
    Nᵢ, SStᵢ, SSeᵢ = anova(Qᵢⱼ...)
    k = length(Nᵢ)
    t3 = VarianceEqualityTest{Chisq}(Nᵢ, SStᵢ, SSeᵢ, k-1, sum(Nᵢ)-k,
            ("Fligner-Killeen test","Variances","FK"))
end
