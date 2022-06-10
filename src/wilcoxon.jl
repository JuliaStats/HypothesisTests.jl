# Wilcoxon.jl
# Wilcoxon signed rank tests
#
# Copyright (C) 2012   Simon Kornblith
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

export SignedRankTest, ExactSignedRankTest, ApproximateSignedRankTest

## COMMON SIGNED RANK

# Automatic exact/normal selection
"""
    SignedRankTest(x::AbstractVector{<:Real})
    SignedRankTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Perform a Wilcoxon signed rank test of the null hypothesis that the distribution of `x`
(or the difference `x - y` if `y` is provided) has zero median against the alternative
hypothesis that the median is non-zero.

When there are no tied ranks and ≤50 samples, or tied ranks and ≤15 samples,
`SignedRankTest` performs an exact signed rank test. In all other cases,
`SignedRankTest` performs an approximate signed rank test. Behavior may be further
controlled by using [`ExactSignedRankTest`](@ref) or [`ApproximateSignedRankTest`](@ref)
directly.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function SignedRankTest(x::AbstractVector{T}) where T<:Real
    (W, ranks, signs, tie_adjustment, n, median) = signedrankstats(x)
    n_nonzero = length(ranks)
    if n_nonzero <= 15 || (n_nonzero <= 50 && tie_adjustment == 0)
        ExactSignedRankTest(x, W, ranks, signs, tie_adjustment, n, median)
    else
        ApproximateSignedRankTest(x, W, ranks, signs, tie_adjustment, n, median)
    end
end
SignedRankTest(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real,S<:Real} = SignedRankTest(x - y)

# Get W and absolute ranks for signed rank test
function signedrankstats(x::AbstractVector{S}) where S<:Real
   nonzero_x = x[x .!= 0]
   (ranks, tieadj) = tiedrank_adj(abs.(nonzero_x))
   W = 0.0
   for i = 1:length(nonzero_x)
       if nonzero_x[i] > 0
           W += ranks[i]
       end
   end
   (W, ranks, nonzero_x .> 0, tieadj, length(x), median(x))
end

## EXACT WILCOXON SIGNED RANK TEST

struct ExactSignedRankTest{T<:Real} <: HypothesisTest
    vals::Vector{T} # original values
    W::Float64              # test statistic: Wilcoxon rank-sum statistic
    ranks::Vector{Float64}           # ranks without ties (zero values)
    signs::BitArray{1}      # signs of input of ranks
    tie_adjustment::Float64 # adjustment for ties
    n::Int                  # number of observations
    median::Float64         # sample median
end
"""
    ExactSignedRankTest(x::AbstractVector{<:Real}[, y::AbstractVector{<:Real}])

Perform a Wilcoxon exact signed rank U test of the null hypothesis that the distribution of
`x` (or the difference `x - y` if `y` is provided) has zero median against the alternative
hypothesis that the median is non-zero.

When there are no tied ranks, the exact p-value is computed using the `psignrank` function
from the `Rmath` package. In the presence of tied ranks, a p-value is computed by exhaustive
enumeration of permutations, which can be very slow for even moderately sized data sets.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
ExactSignedRankTest(x::AbstractVector{T}) where {T<:Real} =
    ExactSignedRankTest(x, signedrankstats(x)...)
ExactSignedRankTest(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real,T<:Real} =
    ExactSignedRankTest(x - y)

testname(::ExactSignedRankTest) = "Exact Wilcoxon signed rank test"
population_param_of_interest(x::ExactSignedRankTest) = ("Location parameter (pseudomedian)", 0, x.median) # parameter of interest: name, value under h0, point estimate
default_tail(test::ExactSignedRankTest) = :both

function show_params(io::IO, x::ExactSignedRankTest, ident)
    println(io, ident, "number of observations:      ", x.n)
    println(io, ident, "Wilcoxon rank-sum statistic: ", x.W)
    print(io, ident, "rank sums:                   ")
    show(io, [sum(x.ranks[x.signs]), sum(x.ranks[map(!, x.signs)])])
    println(io)
    println(io, ident, "adjustment for ties:         ", x.tie_adjustment)
end

# Enumerate all possible Wilcoxon rank-sum results for a given vector, determining left-
# and right-tailed p values
function signedrankenumerate(x::ExactSignedRankTest)
    le = 0
    gr = 0
    n = length(x.ranks)
    tot = 2^n
    for i = 0:tot-1
        # Interpret bits of i as signs to generate wp for all possible sign combinations
        Wp = 0
        b = i
        j = 1
        while b != 0
            Wp += (b & 1)*x.ranks[j]
            j += 1
            b >>= 1
        end
        le += Wp <= x.W
        gr += Wp >= x.W
    end
    (le/tot, gr/tot)
end

function pvalue(x::ExactSignedRankTest; tail=:both)
    check_tail(tail)

    n = length(x.ranks)
    if n == 0
        1.0
    elseif x.tie_adjustment == 0
        # Compute exact p-value using method from Rmath, which is fast but cannot account for ties
        if tail == :both
            if x.W <= n * (n + 1)/4
                2 * psignrank(x.W, n, true)
            else
                2 * psignrank(x.W - 1, n, false)
            end
        elseif tail == :left
            psignrank(x.W, n, true)
        else
            psignrank(x.W - 1, n, false)
        end
    else
        # Compute exact p-value by enumerating all possible ranks in the tied data
        if tail == :both
            min(1, 2 * minimum(signedrankenumerate(x)))
        elseif tail == :left
            first(signedrankenumerate(x))
        else
            last(signedrankenumerate(x))
        end
    end
end

StatsBase.confint(x::ExactSignedRankTest; level::Real=0.95, tail=:both) = calculate_ci(x.vals, level, tail=tail)


## APPROXIMATE SIGNED RANK TEST

struct ApproximateSignedRankTest{T<:Real} <: HypothesisTest
    vals::Vector{T} # original values
    W::Float64              # test statistic: Wilcoxon rank-sum statistic
    ranks::Vector{Float64} # ranks without ties (zero values)
    signs::BitArray{1}      # signs of input of ranks
    tie_adjustment::Float64 # adjustment for ties
    n::Int                  # number of observations
    median::Float64         # sample median
    mu::Float64             # normal approximation: mean
    sigma::Float64          # normal approximation: std
end
"""
    ApproximateSignedRankTest(x::AbstractVector{<:Real}[, y::AbstractVector{<:Real}])

Perform a Wilcoxon approximate signed rank U test of the null hypothesis that the
distribution of `x` (or the difference `x - y` if `y` is provided) has zero median against
the alternative hypothesis that the median is non-zero.

The p-value is computed using a normal approximation to the distribution of the signed rank
statistic:
```math
    \\begin{align*}
        μ & = \\frac{n(n + 1)}{4}\\\\
        σ & = \\frac{n(n + 1)(2 * n + 1)}{24} - \\frac{a}{48}\\\\
        a & = \\sum_{t \\in \\mathcal{T}} t^3 - t
    \\end{align*}
```
where ``\\mathcal{T}`` is the set of the counts of tied values at each tied position.

Implements: [`pvalue`](@ref), [`confint`](@ref)
"""
function ApproximateSignedRankTest(x::Vector, W::Float64, ranks::Vector{T}, signs::BitArray{1}, tie_adjustment::Float64, n::Int, median::Float64) where T<:Real
    nz = length(ranks) # num non-zeros
    mu = W - nz * (nz + 1)/4
    std = sqrt(nz * (nz + 1) * (2 * nz + 1) / 24 - tie_adjustment / 48)
    ApproximateSignedRankTest(x, W, ranks, signs, tie_adjustment, n, median, mu, std)
end
ApproximateSignedRankTest(x::AbstractVector{T}) where {T<:Real} =
    ApproximateSignedRankTest(x, signedrankstats(x)...)
ApproximateSignedRankTest(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real,T<:Real} =
    ApproximateSignedRankTest(x - y)

testname(::ApproximateSignedRankTest) = "Approximate Wilcoxon signed rank test"
population_param_of_interest(x::ApproximateSignedRankTest) = ("Location parameter (pseudomedian)", 0, x.median) # parameter of interest: name, value under h0, point estimate
default_tail(test::ApproximateSignedRankTest) = :both

function show_params(io::IO, x::ApproximateSignedRankTest, ident)
    println(io, ident, "number of observations:      ", x.n)
    println(io, ident, "Wilcoxon rank-sum statistic: ", x.W)
    print(io, ident, "rank sums:                   ")
    show(io, [sum(x.ranks[x.signs]), sum(x.ranks[map(!, x.signs)])])
    println(io)
    println(io, ident, "adjustment for ties:         ", x.tie_adjustment)
    println(io, ident, "normal approximation (μ, σ): ", (x.mu, x.sigma))
end

function pvalue(x::ApproximateSignedRankTest; tail=:both)
    check_tail(tail)

    if x.mu == x.sigma == 0
        1.0
    elseif tail == :both
        2 * ccdf(Normal(), abs(x.mu - 0.5 * sign(x.mu))/x.sigma)
    elseif tail == :left
        cdf(Normal(), (x.mu + 0.5)/x.sigma)
    else # tail == :right
        ccdf(Normal(), (x.mu - 0.5)/x.sigma)
    end
end

StatsBase.confint(x::ApproximateSignedRankTest; level::Real=0.95, tail=:both) = calculate_ci(x.vals, level, tail=tail)

# implementation method inspired by these notes: http://www.stat.umn.edu/geyer/old03/5102/notes/rank.pdf
function calculate_ci(x::AbstractVector, level::Real=0.95; tail=:both)
    check_level(level)
    check_tail(tail)

    if tail == :both
        c = level
    else
        c = 1 - 2 * (1-level)
    end
    n = length(x)
    m = div(n * (n + 1), 2)
    k_range = 1:div(m, 2)
    l = [1 - 2 * psignrank(i, n, true) for i in k_range]
    k = argmin(abs.(l .- c))
    vals = Float64[]
    enumerated = enumerate(x)
    for (outer_index, outer_value) in enumerated
        for (inner_index, inner_value) in enumerated
            if outer_index > inner_index
                continue
            end
            push!(vals, (inner_value + outer_value) / 2)
        end
    end
    sort!(vals)
    left = vals[k + 1]
    right = vals[m - k]
    if tail == :both
        return (left, right)
    elseif tail == :left
        return (left, Inf)
    else # tail == :right
        return (-Inf, right)
    end
end
