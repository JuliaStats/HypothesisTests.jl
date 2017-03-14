# Wilcoxon.jl
# Wilcoxon rank sum (Mann-Whitney U) and signed rank tests in Julia
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
function SignedRankTest{T<:Real}(x::AbstractVector{T})
    (W, ranks, signs, tie_adjustment, n, median) = signedrankstats(x)
    n_nonzero = length(ranks)
    if n_nonzero <= 15 || (n_nonzero <= 50 && tie_adjustment == 0)
        ExactSignedRankTest(x, W, ranks, signs, tie_adjustment, n, median)
    else
        ApproximateSignedRankTest(x, W, ranks, signs, tie_adjustment, n, median)
    end
end
SignedRankTest{T<:Real,S<:Real}(x::AbstractVector{T}, y::AbstractVector{S}) = SignedRankTest(x - y)

# Get W and absolute ranks for signed rank test
function signedrankstats{S<:Real}(x::AbstractVector{S})
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

immutable ExactSignedRankTest{T<:Real} <: HypothesisTest
    vals::Vector{T} # original values
    W::Float64              # test statistic: Wilcoxon rank-sum statistic
    ranks::Vector{Float64}           # ranks without ties (zero values)
    signs::BitArray{1}      # signs of input of ranks
    tie_adjustment::Float64 # adjustment for ties
    n::Int                  # number of observations
    median::Float64         # sample median
end
ExactSignedRankTest{T<:Real}(x::AbstractVector{T}) =
    ExactSignedRankTest(x, signedrankstats(x)...)
ExactSignedRankTest{S<:Real,T<:Real}(x::AbstractVector{S}, y::AbstractVector{T}) =
    ExactSignedRankTest(x - y)

testname(::ExactSignedRankTest) = "Exact Wilcoxon signed rank test"
population_param_of_interest(x::ExactSignedRankTest) = ("Location parameter (pseudomedian)", 0, x.median) # parameter of interest: name, value under h0, point estimate

function show_params(io::IO, x::ExactSignedRankTest, ident)
    println(io, ident, "number of observations:      ", x.n)
    println(io, ident, "Wilcoxon rank-sum statistic: ", x.W)
    println(io, ident, "rank sums:                   ", [sum(x.ranks[x.signs]), sum(x.ranks[map(!, x.signs)])])
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
    n = length(x.ranks)
    if n == 0
        1
    elseif x.tie_adjustment == 0
        # Compute exact p-value using method from Rmath, which is fast but cannot account for ties
        if tail == :both
            if x.W <= n * (n + 1)/4
                p = 2 * psignrank(x.W, n, true)
            else
                p = 2 * psignrank(x.W - 1, n, false)
            end
        elseif tail == :left
            psignrank(x.W, n, true)
        elseif tail == :right
            psignrank(x.W - 1, n, false)
        end
    else
        # Compute exact p-value by enumerating all possible ranks in the tied data
        if tail == :both
            min(1, 2 * minimum(signedrankenumerate(x)))
        elseif tail == :left
            singedrankenumerate(x.W, x.ranks)[1]
        elseif tail == :right
            signedrankenumerate(x.W, x.ranks)[2]
        end
    end
end

StatsBase.confint(x::ExactSignedRankTest, alpha::Real=0.05; tail=:both) = calculate_ci(x.vals, alpha; tail=tail)


## APPROXIMATE SIGNED RANK TEST

immutable ApproximateSignedRankTest{T<:Real} <: HypothesisTest
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

function ApproximateSignedRankTest{T<:Real}(x::Vector, W::Float64, ranks::Vector{T}, signs::BitArray{1}, tie_adjustment::Float64, n::Int, median::Float64)
    nz = length(ranks) # num non-zeros
    mu = W - nz * (nz + 1)/4
    std = sqrt(nz * (nz + 1) * (2 * nz + 1) / 24 - tie_adjustment / 48)
    ApproximateSignedRankTest(x, W, ranks, signs, tie_adjustment, n, median, mu, std)
end
ApproximateSignedRankTest{T<:Real}(x::AbstractVector{T}) =
    ApproximateSignedRankTest(x, signedrankstats(x)...)
ApproximateSignedRankTest{S<:Real,T<:Real}(x::AbstractVector{S}, y::AbstractVector{T}) =
    ApproximateSignedRankTest(x - y)

testname(::ApproximateSignedRankTest) = "Approximate Wilcoxon signed rank test"
population_param_of_interest(x::ApproximateSignedRankTest) = ("Location parameter (pseudomedian)", 0, x.median) # parameter of interest: name, value under h0, point estimate

function show_params(io::IO, x::ApproximateSignedRankTest, ident)
    println(io, ident, "number of observations:      ", x.n)
    println(io, ident, "Wilcoxon rank-sum statistic: ", x.W)
    println(io, ident, "rank sums:                   ", [sum(x.ranks[x.signs]), sum(x.ranks[map(!, x.signs)])])
    println(io, ident, "adjustment for ties:         ", x.tie_adjustment)
    println(io, ident, "normal approximation (μ, σ): ", (x.mu, x.sigma))
end

function pvalue(x::ApproximateSignedRankTest; tail=:both)
    if x.mu == x.sigma == 0
        1
    else
        if tail == :both
            2 * ccdf(Normal(), abs(x.mu - 0.5 * sign(x.mu))/x.sigma)
        elseif tail == :left
            cdf(Normal(), (x.mu + 0.5)/x.sigma)
        elseif tail == :right
            ccdf(Normal(), (x.mu - 0.5)/x.sigma)
        end
    end
end

StatsBase.confint(x::ApproximateSignedRankTest, alpha::Real=0.05; tail=:both) = calculate_ci(x.vals, alpha; tail=tail)

# implementation method inspired by these notes: http://www.stat.umn.edu/geyer/old03/5102/notes/rank.pdf
function calculate_ci(x::AbstractVector, alpha::Real=0.05; tail=:both)
    check_alpha(alpha)

    if tail == :both
        c = 1 - alpha
    else
        c = 1 - 2 * alpha
    end
    n = length(x)
    m = div(n * (n + 1), 2)
    k_range = 1:div(m, 2)
    l = [1 - 2 * psignrank(i, n, true) for i in k_range]
    k = indmin(abs.(l-c))
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
    elseif tail == :right
        return (-Inf, right)
    end
end
