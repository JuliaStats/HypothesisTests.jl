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
    SignedRankTest(x::AbstractVector{<:Real}; method = :auto)
    SignedRankTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; method = :auto)

Perform a Wilcoxon signed rank test of the null hypothesis that the distribution of `x`
(or the difference `x - y` if `y` is provided) is symmetric about zero, against the
alternative that it is not. Under a location model, where the distribution is symmetric about
some `θ`, that is `θ = 0` against `θ ≠ 0`, with `tail = :left` and `tail = :right` giving
`θ < 0` and `θ > 0`.

Symmetry is what the test needs, not a zero median: against a null of zero median alone it
does not hold its level. The location it estimates is the pseudomedian (see
[`hodgeslehmann`](@ref)), which equals the median when the distribution is symmetric.

`SignedRankTest` chooses between the exact and the approximate test by the tie pattern and
the number `n` of non-zero observations: with no tied ranks it is exact for `n ≤ 50`, with
tied ranks for `n ≤ 15`, and approximate above whichever of those two applies.

`method` overrides that choice:

  - `:auto` (the default) applies the rule above.
  - `:exact` and `:approximate` force the corresponding test, which is what an analysis
    that must reproduce across versions of this package should do.
  - a callable is passed `(; n, n_nonzero, ties, tie_adjustment)` and must return
    `:exact` or `:approximate`.

Equivalently, construct [`ExactSignedRankTest`](@ref) or
[`ApproximateSignedRankTest`](@ref) directly.

Implements: [`pvalue`](@ref), [`confint`](@ref), [`hodgeslehmann`](@ref)
"""
function SignedRankTest(x::AbstractVector{T}; method = :auto) where T<:Real
    (W, ranks, signs, tie_adjustment, n, median) = signedrankstats(x)
    # `ranks` covers the non-zero observations alone, so its length is that count
    n_nonzero = length(ranks)
    # the named tuple the `method` callable is documented to receive, and the automatic
    # rule it defaults to, which is the threshold this constructor has always applied
    stats = (n = n, n_nonzero = n_nonzero, ties = tie_adjustment != 0,
             tie_adjustment = tie_adjustment)
    default = n_nonzero <= 15 || (n_nonzero <= 50 && tie_adjustment == 0) ? :exact : :approximate
    if resolve_rank_method(method, stats, default) === :exact
        ExactSignedRankTest(x, W, ranks, signs, tie_adjustment, n, median)
    else
        ApproximateSignedRankTest(x, W, ranks, signs, tie_adjustment, n, median)
    end
end
SignedRankTest(x::AbstractVector{T}, y::AbstractVector{S}; method = :auto) where {T<:Real,S<:Real} =
    SignedRankTest(x - y; method = method)

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

struct ExactSignedRankTest <: HypothesisTest
    vals::Vector{Float64}   # original values
    W::Float64              # test statistic: the signed rank statistic W+
    ranks::Vector{Float64}           # midranks of |d| over the non-zero observations
    signs::BitArray{1}      # signs of input of ranks
    tie_adjustment::Float64 # adjustment for ties
    n::Int                  # number of observations
    median::Float64         # sample median
end
"""
    ExactSignedRankTest(x::AbstractVector{<:Real}[, y::AbstractVector{<:Real}])

Perform an exact Wilcoxon signed rank test of the null hypothesis that the distribution of
`x` (or the difference `x - y` if `y` is provided) is symmetric about zero, against the
alternative that it is not. See [`SignedRankTest`](@ref) on what that null does and does not
assume, and on the alternatives the `tail` keyword selects.

When there are no tied ranks, the exact p-value is computed using the `signrankcdf` and `signrankccdf`
functions from the `StatsFuns` package. In the presence of tied ranks, a p-value is computed by exhaustive
enumeration of the ``2^n`` sign assignments over the non-zero observations.

The tied route is bounded by [`MAX_EXACT_ENUMERATION_N`](@ref): beyond it this test
refuses rather than enumerate indefinitely, and `method = :approximate` is the way on.

Implements: [`pvalue`](@ref), [`confint`](@ref), [`hodgeslehmann`](@ref)
"""
ExactSignedRankTest(x::AbstractVector{T}) where {T<:Real} =
    ExactSignedRankTest(x, signedrankstats(x)...)
ExactSignedRankTest(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real,T<:Real} =
    ExactSignedRankTest(x - y)

testname(::ExactSignedRankTest) = "Exact Wilcoxon signed rank test"
population_param_of_interest(x::ExactSignedRankTest) = ("Location parameter (pseudomedian)", 0, hodgeslehmann(x)) # parameter of interest: name, value under h0, point estimate
default_tail(test::ExactSignedRankTest) = :both

function show_params(io::IO, x::ExactSignedRankTest, ident)
    println(io, ident, "number of observations:      ", x.n)
    println(io, ident, "non-zero observations:       ", length(x.ranks))
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
    check_exact_enumeration(n)
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

function StatsAPI.pvalue(x::ExactSignedRankTest; tail=:both)
    check_tail(tail)

    n = length(x.ranks)
    if n == 0
        1.0
    elseif x.tie_adjustment == 0
        # Compute exact p-value using method from StatsFuns, which is fast but cannot account for ties
        if tail == :both
            # The smaller tail, doubled and clipped. Doubling a discrete tail can
            # overshoot: where n(n+1)/2 is even the null mean is attainable, and at
            # W == n(n+1)/4 each tail is (1 + P(W == mean))/2, so the doubling gives
            # 1 + P(W == mean), which is 1.25 for n = 3 at W = 3. The tied branch below
            # and both exact Mann-Whitney branches clip for the same reason.
            p = x.W <= n * (n + 1)/4 ? signrankcdf(n, x.W) : signrankccdf(n, x.W - 1)
            min(2 * p, 1.0)
        elseif tail == :left
            signrankcdf(n, x.W)
        else
            signrankccdf(n, x.W - 1)
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

# The Walsh averages that the interval and the point estimate are read off.
#
# `signedrankstats` drops zero differences before ranking, so the statistic and the
# p-value describe the non-zero differences; the interval and the estimate must
# describe the same sample. R's `wilcox.test` drops them likewise.
function signedrank_pairwise_estimates(vals::AbstractVector{<:Real})
    nonzero = filter(!iszero, vals)
    # every difference is zero: nothing to drop, and the pairwise estimates degenerate to
    # the point zero either way
    return walsh_averages(isempty(nonzero) ? vals : nonzero)
end

hodgeslehmann(x::ExactSignedRankTest) = median(signedrank_pairwise_estimates(x.vals))

# Ties are a caveat here rather than a correction: under ties the null distribution of
# the statistic is no longer the untied one this inverts, so the achieved coverage is
# approximate. (R declines to compute an exact interval at all in that case and falls
# back to the normal approximation, which on the samples tested lands on the same order
# statistics.) `ApproximateSignedRankTest` does account for ties, through its variance.
function StatsAPI.confint(x::ExactSignedRankTest; level::Real=0.95, tail=:both)
    alpha = ci_alpha(level, tail)
    n = length(x.ranks)
    vals = signedrank_pairwise_estimates(x.vals)
    # this route still inverts the exact distribution, a lattice recursion per candidate
    # endpoint, so it keeps the bound the scan it replaces was given
    check_exact_ci_cost(length(vals),
        "Pass `method = :approximate` for the normal-approximation interval, which " *
        "inverts a closed form and is bounded by memory alone.")
    k = exact_ci_index(length(vals), alpha, i -> signrankcdf(n, i))
    return ci_from_estimates(vals, k, tail)
end


## APPROXIMATE SIGNED RANK TEST

struct ApproximateSignedRankTest <: HypothesisTest
    vals::Vector{Float64}   # original values
    W::Float64              # test statistic: the signed rank statistic W+
    ranks::Vector{Float64} # midranks of |d| over the non-zero observations
    signs::BitArray{1}      # signs of input of ranks
    tie_adjustment::Float64 # adjustment for ties
    n::Int                  # number of observations
    median::Float64         # sample median
    mu::Float64             # normal approximation: mean
    sigma::Float64          # normal approximation: std
end
"""
    ApproximateSignedRankTest(x::AbstractVector{<:Real}[, y::AbstractVector{<:Real}])

Perform an approximate Wilcoxon signed rank test of the null hypothesis that the
distribution of `x` (or the difference `x - y` if `y` is provided) is symmetric about zero,
against the alternative that it is not. See
[`SignedRankTest`](@ref) on what that null does and does not assume.

The p-value is computed using a normal approximation to the distribution of the signed rank
statistic ``W^+``, which under the null has mean and variance
```math
    \\begin{align*}
        μ_0 & = \\frac{n(n + 1)}{4}\\\\
        σ^2 & = \\frac{n(n + 1)(2n + 1)}{24} - \\frac{a}{48}\\\\
        a & = \\sum_{t \\in \\mathcal{T}} t^3 - t
    \\end{align*}
```
where ``\\mathcal{T}`` is the set of the counts of tied values at each tied position and
``n`` counts the non-zero observations. What `show` reports as `normal approximation (μ, σ)`
is the pair ``(W^+ - μ_0, σ)``: the statistic centred at its null mean, and the
tie-corrected standard deviation, not ``μ_0`` itself.

The confidence interval inverts the same approximation, rather than the exact null
distribution.

Implements: [`pvalue`](@ref), [`confint`](@ref), [`hodgeslehmann`](@ref)
"""
function ApproximateSignedRankTest(x::Vector, W::Float64, ranks::Vector{T}, signs::BitArray{1}, tie_adjustment::Float64, n::Int, median::Real) where T<:Real
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
population_param_of_interest(x::ApproximateSignedRankTest) = ("Location parameter (pseudomedian)", 0, hodgeslehmann(x)) # parameter of interest: name, value under h0, point estimate
default_tail(test::ApproximateSignedRankTest) = :both

function show_params(io::IO, x::ApproximateSignedRankTest, ident)
    println(io, ident, "number of observations:      ", x.n)
    println(io, ident, "non-zero observations:       ", length(x.ranks))
    println(io, ident, "Wilcoxon rank-sum statistic: ", x.W)
    print(io, ident, "rank sums:                   ")
    show(io, [sum(x.ranks[x.signs]), sum(x.ranks[map(!, x.signs)])])
    println(io)
    println(io, ident, "adjustment for ties:         ", x.tie_adjustment)
    println(io, ident, "normal approximation (μ, σ): ", (x.mu, x.sigma))
end

function StatsAPI.pvalue(x::ApproximateSignedRankTest; tail=:both)
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

hodgeslehmann(x::ApproximateSignedRankTest) = median(signedrank_pairwise_estimates(x.vals))

# The exact null distribution is not consulted here: an approximate test gets an
# approximate interval, from the same normal approximation (mean, and variance
# corrected for ties) that its p-value uses.
function StatsAPI.confint(x::ApproximateSignedRankTest; level::Real=0.95, tail=:both)
    alpha = ci_alpha(level, tail)
    vals = signedrank_pairwise_estimates(x.vals)
    m = length(vals)
    k = normal_ci_index(m, m / 2, x.sigma, alpha)
    return ci_from_estimates(vals, k, tail)
end
