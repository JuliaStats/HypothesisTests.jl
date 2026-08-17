# Wilcoxon.jl
# Wilcoxon rank sum (Mann-Whitney U) tests
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

export MannWhitneyUTest, ExactMannWhitneyUTest, ApproximateMannWhitneyUTest

## COMMON MANN-WHITNEY U

# Automatic exact/normal selection
"""
    MannWhitneyUTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; method = :auto)

Perform a Mann-Whitney U test of the null hypothesis that `x` and `y`
are drawn from the same distribution, against the alternative that one tends to exceed
the other.

Equality of the two distributions is what makes the test exact: it renders the pooled
observations exchangeable, which is what the null distribution of the statistic rests on.
This is *not* a test of `P(x > y) = P(y > x)` alone. That weaker equality permits the two
distributions to differ in spread, and then the test does not hold its level: for
`Normal(0, 1)` against `Normal(0, 9)`, where it holds exactly, a nominal 0.05 test rejects
about 13% of the time at `nx, ny = 30, 10`, and about 1.6% at `10, 30`.

Under a shift model the location estimated is the median of `x - y`, reported by
[`hodgeslehmann`](@ref).

The Mann-Whitney U test is sometimes known as the Wilcoxon rank-sum test.

`MannWhitneyUTest` chooses between the exact and the approximate test by the tie pattern and
the pooled sample size `nx + ny`: with no tied ranks it is exact for `nx + ny ≤ 50`, with tied
ranks for `nx + ny ≤ 10`, and approximate above whichever of those two applies.

`method` overrides that choice:

  - `:auto` (the default) applies the rule above.
  - `:exact` and `:approximate` force the corresponding test, which is what an analysis
    that must reproduce across versions of this package should do.
  - a callable is passed `(; nx, ny, ties, tie_adjustment)` and must return `:exact` or
    `:approximate`.

Equivalently, construct [`ExactMannWhitneyUTest`](@ref) or
[`ApproximateMannWhitneyUTest`](@ref) directly.

Implements: [`pvalue`](@ref), [`confint`](@ref), [`hodgeslehmann`](@ref)
"""
function MannWhitneyUTest(x::AbstractVector{S}, y::AbstractVector{T}; method = :auto) where {S<:Real,T<:Real}
    fields = mwustats(x, y)
    (U, ranks, tieadj, nx, ny, median) = fields
    # the named tuple the `method` callable is documented to receive, and the automatic
    # rule it defaults to, which is the threshold this constructor has always applied
    stats = (nx = nx, ny = ny, ties = tieadj != 0, tie_adjustment = tieadj)
    default = nx + ny <= 10 || (nx + ny <= 50 && tieadj == 0) ? :exact : :approximate
    if resolve_rank_method(method, stats, default) === :exact
        ExactMannWhitneyUTest(fields...)
    else
        ApproximateMannWhitneyUTest(fields...)
    end
end

# Get U, ranks, and tie adjustment for Mann-Whitney U test
# U is the sum of adjusted ranks in the first sample minus the minimal sum of ranks (ie sum(1:length(x))
# The samples themselves are carried through as well: `confint` and `hodgeslehmann`
# need the cross-group differences, which cannot be recovered from the ranks. They are
# returned in the caller's own element type and narrowed to `Float64` by the struct
# fields, so ranking and the median are computed at the input's precision and only what
# is stored is converted.
function mwustats(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    ranks, tieadj = tiedrank_adj([x; y])
    nx = length(x)
    ny = length(y)
    if nx <= ny
        U = sum(@view(ranks[begin:(begin + (nx-1))])) - nx*(nx+1)/2
    else
        # Sum of adjusted ranks of first and second sample sums to (nx + ny)*(nx + ny + 1)/2, hence
        # U = (nx + ny)*(nx + ny + 1)/2 - sum(ranks_y) - nx*(nx + 1)/2
        U = ny*(2 * nx + ny + 1)/2 - sum(@view(ranks[(begin + nx):end]))
    end
    return (U, ranks, tieadj, nx, ny, median(x)-median(y), x, y)
end


## EXACT MANN-WHITNEY U TEST
struct ExactMannWhitneyUTest <: HypothesisTest
    U::Float64              # test statistic: Mann-Whitney-U statistic
    ranks::Vector{Float64}  # ranks
    tie_adjustment::Float64 # adjustment for ties
    nx::Int                 # number of observations
    ny::Int
    median::Float64         # difference of sample medians
    x::Vector{Float64}      # original values, first sample
    y::Vector{Float64}      # original values, second sample
end

"""
    ExactMannWhitneyUTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Perform an exact Mann-Whitney U test of the null hypothesis that `x` and `y` are drawn
from the same distribution, against the alternative that one tends to exceed the other.
See [`MannWhitneyUTest`](@ref) on what that null does and does not assume.

When there are no tied ranks, the exact p-value is computed using the `wilcoxcdf` and `wilcoxccdf`
functions from the `StatsFuns` package. In the presence of tied ranks, a p-value is computed by exhaustive
enumeration of the assignments of the pooled midranks to the smaller sample.

The tied route is bounded by [`MAX_EXACT_ENUMERATION_N`](@ref): beyond it this test
refuses rather than enumerate indefinitely, and `method = :approximate` is the way on.

Implements: [`pvalue`](@ref), [`confint`](@ref), [`hodgeslehmann`](@ref)
"""
ExactMannWhitneyUTest(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real,T<:Real} =
    ExactMannWhitneyUTest(mwustats(x, y)...)

testname(::ExactMannWhitneyUTest) = "Exact Mann-Whitney U test"
population_param_of_interest(x::ExactMannWhitneyUTest) = ("Location shift", 0, hodgeslehmann(x)) # parameter of interest: name, value under h0, point estimate
default_tail(test::ExactMannWhitneyUTest) = :both

function show_params(io::IO, x::ExactMannWhitneyUTest, ident)
    print(io, ident, "number of observations in each group: ")
    show(io, [x.nx, x.ny])
    println(io)
    println(io, ident, "Mann-Whitney-U statistic:             ", x.U)
    print(io, ident, "rank sums:                            ")
    show(io, [sum(@view(x.ranks[begin:(begin + (x.nx - 1))])), sum(@view(x.ranks[(begin + x.nx):end]))])
    println(io)
    println(io, ident, "adjustment for ties:                  ", x.tie_adjustment)
end

# Enumerate all possible Mann-Whitney U results for a given vector,
# determining left-and right-tailed p values
function mwuenumerate(x::ExactMannWhitneyUTest)
    # Compute sum of ranks of the smaller group
    nx = x.nx
    ny = x.ny
    check_exact_enumeration(nx, ny)
    if nx <= ny
        n = nx
        R = x.U + n*(n + 1)/2
    else
        n = ny
        R = ny*(2 * nx + ny + 1)/2 - x.U
    end

    # Compare with sum of ranks of all possible groups of the same size
    le = gr = tot = 0
    for comb in combinations(x.ranks, n)
        Rp = sum(comb)
        tot += 1
        le += Rp <= R
        gr += Rp >= R
    end

    # Adjust "less"/"greater" in case the first group was not the smaller one
    if nx > n
        le, gr = gr, le
    end

    (le/tot, gr/tot)
end

function StatsAPI.pvalue(x::ExactMannWhitneyUTest; tail=:both)
    check_tail(tail)

    if x.tie_adjustment == 0
        # Compute exact p-value using method from StatsFuns, which is fast but
        # cannot account for ties
        if tail == :both
            p = wilcoxcdf(x.nx, x.ny, min(x.U, x.nx * x.ny - x.U))
            min(2 * p, 1.0)
        elseif tail == :left
            wilcoxcdf(x.nx, x.ny, x.U)
        else # tail == :right
            wilcoxccdf(x.nx, x.ny, x.U - 1)
        end
    else
        # Compute exact p-value by enumerating possible ranks in the tied data
        if tail == :both
            min(1.0, 2 * minimum(mwuenumerate(x)))
        elseif tail == :left
            mwuenumerate(x)[1]
        else # tail == :right
            mwuenumerate(x)[2]
        end
    end
end

hodgeslehmann(x::ExactMannWhitneyUTest) = median(cross_differences(x.x, x.y))

# Under ties the achieved coverage is approximate rather than exact: the null
# distribution inverted here is the untied one, where the tied data has the
# conditional distribution `mwuenumerate` computes. The classical construction
# inverts the untied one anyway; R declines an exact interval in that case.
function StatsAPI.confint(x::ExactMannWhitneyUTest; level::Real=0.95, tail=:both)
    alpha = ci_alpha(level, tail)
    # before the differences are formed, so an oversized request is refused rather than
    # paid for: `exact_ci_index` below runs a lattice recursion per candidate endpoint
    check_exact_ci_cost(x.nx * x.ny,
        "Pass `method = :approximate` for the normal-approximation interval, which " *
        "inverts a closed form and is bounded by memory alone.")
    vals = cross_differences(x.x, x.y)
    k = exact_ci_index(length(vals), alpha, u -> wilcoxcdf(x.nx, x.ny, u); tail = tail)
    return ci_from_estimates(vals, k, tail)
end

struct ApproximateMannWhitneyUTest <: HypothesisTest
    U::Float64              # test statistic: Mann-Whitney-U statistic
    ranks::Vector{Float64}  # ranks
    tie_adjustment::Float64 # adjustment for ties
    nx::Int                 # number of observations
    ny::Int
    median::Float64         # difference of sample medians
    mu::Float64             # normal approximation: mean
    sigma::Float64          # normal approximation: std
    x::Vector{Float64}      # original values, first sample
    y::Vector{Float64}      # original values, second sample
end

## APPROXIMATE MANN-WHITNEY U TEST
"""
    ApproximateMannWhitneyUTest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Perform an approximate Mann-Whitney U test of the null hypothesis that `x` and `y` are drawn
from the same distribution, against the alternative that one tends to exceed the other.
See [`MannWhitneyUTest`](@ref) on what that null does and does not assume.

The p-value is computed using a normal approximation to the distribution of the
Mann-Whitney U statistic, which under the null has mean and variance
```math
    \\begin{align*}
        μ_0 & = \\frac{n_x n_y}{2}\\\\
        σ^2 & = \\frac{n_x n_y}{12}\\left(n_x + n_y + 1 - \\frac{a}{(n_x + n_y)(n_x +
            n_y - 1)}\\right)\\\\
        a & = \\sum_{t \\in \\mathcal{T}} t^3 - t
    \\end{align*}
```
where ``\\mathcal{T}`` is the set of the counts of tied values at each tied position.
What `show` reports as `normal approximation (μ, σ)` is the pair ``(U - μ_0, σ)``: the
statistic centred at its null mean, and the tie-corrected standard deviation, not ``μ_0``
itself.

The confidence interval inverts the same approximation, rather than the exact null
distribution.

Implements: [`pvalue`](@ref), [`confint`](@ref), [`hodgeslehmann`](@ref)
"""
function ApproximateMannWhitneyUTest(U::Real, ranks::AbstractVector{<:Real},
    tie_adjustment::Real, nx::Int, ny::Int, median::Real,
    x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = nx + ny
    nxny = nx * ny
    mu = U - nxny / 2
    sigma = sqrt(nxny / 12 * (n + 1 - tie_adjustment / (n * (n - 1))))
    ApproximateMannWhitneyUTest(U, ranks, tie_adjustment, nx, ny, median, mu, sigma, x, y)
end
ApproximateMannWhitneyUTest(x::AbstractVector{S}, y::AbstractVector{T}) where {S<:Real,T<:Real} =
    ApproximateMannWhitneyUTest(mwustats(x, y)...)

testname(::ApproximateMannWhitneyUTest) = "Approximate Mann-Whitney U test"
population_param_of_interest(x::ApproximateMannWhitneyUTest) = ("Location shift", 0, hodgeslehmann(x)) # parameter of interest: name, value under h0, point estimate
default_tail(test::ApproximateMannWhitneyUTest) = :both

function show_params(io::IO, x::ApproximateMannWhitneyUTest, ident)
    println(io, ident, "number of observations in each group: ", [x.nx, x.ny])
    println(io, ident, "Mann-Whitney-U statistic:             ", x.U)
    println(io, ident, "rank sums:                            ", [sum(@view(x.ranks[begin:(begin + (x.nx - 1))])), sum(@view(x.ranks[(begin + x.nx):end]))])
    println(io, ident, "adjustment for ties:                  ", x.tie_adjustment)
    println(io, ident, "normal approximation (μ, σ):          ", (x.mu, x.sigma))
end

function StatsAPI.pvalue(x::ApproximateMannWhitneyUTest; tail=:both)
    check_tail(tail)

    if x.mu == x.sigma == 0
        1.0
    elseif tail == :both
        p = 2 * ccdf(Normal(), abs(x.mu - 0.5 * sign(x.mu))/x.sigma)
        min(p, 1.0)
    elseif tail == :left
        cdf(Normal(), (x.mu + 0.5)/x.sigma)
    else # tail == :right
        ccdf(Normal(), (x.mu - 0.5)/x.sigma)
    end
end

hodgeslehmann(x::ApproximateMannWhitneyUTest) = median(cross_differences(x.x, x.y))

# `x.mu` is the centred statistic, not the null mean; the null mean of U is nx*ny/2.
function StatsAPI.confint(x::ApproximateMannWhitneyUTest; level::Real=0.95, tail=:both)
    alpha = ci_alpha(level, tail)
    vals = cross_differences(x.x, x.y)
    m = length(vals)
    k = normal_ci_index(m, m / 2, x.sigma, alpha; tail = tail)
    return ci_from_estimates(vals, k, tail)
end
