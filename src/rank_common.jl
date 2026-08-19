# rank_common.jl
# Machinery shared by the Wilcoxon signed rank and Mann-Whitney U tests:
# selection of the exact/approximate method, the pairwise estimates that the
# distribution-free intervals and the point estimates are read off, the two index
# rules used to invert the test, and bounds on the cost of both.

export hodgeslehmann

## METHOD SELECTION

"""
    resolve_rank_method(method, stats::NamedTuple, default::Symbol) -> Symbol

Resolve the `method` keyword of a rank test to either `:exact` or `:approximate`.

`method` may be `:auto` (use `default`, which is the automatic rule of the test in
question), `:exact`, `:approximate`, or a callable mapping `stats` to one of the
latter two.
"""
function resolve_rank_method(method, stats::NamedTuple, default::Symbol)
    if method === :auto
        return default
    elseif method === :exact || method === :approximate
        return method
    elseif method isa Function
        resolved = method(stats)
        if resolved !== :exact && resolved !== :approximate
            throw(ArgumentError(
                "a `method` callable must return `:exact` or `:approximate`, got $(repr(resolved))"))
        end
        return resolved::Symbol
    else
        throw(ArgumentError(
            "`method` must be `:auto`, `:exact`, `:approximate`, or a callable, got $(repr(method))"))
    end
end

## GUARDS ON THE EXACT ROUTES
#
# The exact null distributions themselves need no guard: `StatsFuns.signrankcdf`
# and `StatsFuns.wilcoxcdf` accumulated lattice counts in `Int` and were silently
# invalid once those overflowed (from n = 72, and once binomial(nx + ny, nx)
# exceeded about 3.0e18), but JuliaStats/StatsFuns.jl#221 folded the normaliser
# into the recursions, so the counts are never formed. The `[compat]` floor of
# StatsFuns 2.2.1 is what makes that a settled question rather than a bound here.
# See JuliaStats/StatsFuns.jl#219 and #220.
#
# What does still need bounding is this package's own tied-data enumeration, and
# the pairwise estimates the intervals are read off. Both are cost, not correctness.

"""
Bound on the tied-data enumerations, which is where an exact rank test computes its p-value
when the tie total is not zero.

The signed rank test enumerates ``2^n`` sign assignments over its ``n`` non-zero
observations, and refuses past `n = MAX_EXACT_ENUMERATION_N`. The Mann-Whitney test
enumerates ``\\binom{n_x + n_y}{\\min(n_x, n_y)}`` rank assignments, so the bound is on that
count rather than on the sample size: it refuses once there are more than
`2^MAX_EXACT_ENUMERATION_N` of them. Either way the p-value raises rather than run for an
unbounded time, and `method = :approximate` is the way on; `confint` is unaffected, since it
inverts the tie-free distribution whichever route the p-value took. See #7.
"""
const MAX_EXACT_ENUMERATION_N = 25

function check_exact_enumeration(n::Integer)
    n <= MAX_EXACT_ENUMERATION_N && return nothing
    throw(ArgumentError(
        "an exact test with ties enumerates all 2^n sign assignments, which is not " *
        "feasible for n = $n non-zero observations (limit " *
        "MAX_EXACT_ENUMERATION_N = $MAX_EXACT_ENUMERATION_N). Pass `method = :approximate`."))
end

function check_exact_enumeration(nx::Integer, ny::Integer)
    count = rank_binomial(nx, ny)
    (count !== nothing && count <= Int128(2)^MAX_EXACT_ENUMERATION_N) && return nothing
    throw(ArgumentError(
        "an exact test with ties enumerates all binomial(nx + ny, min(nx, ny)) rank " *
        "assignments, which is not feasible for nx = $nx, ny = $ny (limit " *
        "2^MAX_EXACT_ENUMERATION_N = $(Int128(2)^MAX_EXACT_ENUMERATION_N)). " *
        "Pass `method = :approximate`."))
end

# `binomial(nx + ny, nx)`, or `nothing` when it is too large to represent, which is
# itself an answer, since every bound we compare it against is far smaller.
function rank_binomial(nx::Integer, ny::Integer)
    try
        return binomial(Int128(nx + ny), Int128(min(nx, ny)))
    catch err
        err isa OverflowError && return nothing
        rethrow()
    end
end

## PAIRWISE ESTIMATES

"""
Largest set of pairwise estimates (Walsh averages, or cross-group differences) that will be
materialised. The estimators and intervals below are ``O(n^2)`` in memory; beyond
this bound they refuse rather than exhaust the machine. Displaying such a test still
works: `show` drops its confidence interval line rather than failing to print at all.
See #7 and #97.
"""
const MAX_PAIRWISE_ESTIMATES = 100_000_000

function check_estimate_count(m::Integer, what::AbstractString)
    m <= MAX_PAIRWISE_ESTIMATES && return nothing
    throw(ArgumentError(
        "refusing to form $m $what: more than MAX_PAIRWISE_ESTIMATES = $MAX_PAIRWISE_ESTIMATES. " *
        "The distribution-free interval and the Hodges-Lehmann estimator are computed " *
        "by materialising this set, which is quadratic in the sample size (see #7)."))
end

"""
    walsh_averages(x) -> Vector

The ``n(n+1)/2`` Walsh averages ``(xᵢ + xⱼ)/2`` for ``i ≤ j``, sorted ascending.
Their median is the one-sample Hodges-Lehmann estimator, and the distribution-free
interval for the pseudomedian is a pair of their order statistics.
"""
function walsh_averages(x::AbstractVector{T}) where T<:Real
    n = length(x)
    m = div(n * (n + 1), 2)
    check_estimate_count(m, "Walsh averages")
    out = Vector{float(T)}(undef, m)
    k = 0
    @inbounds for i in 1:n, j in i:n
        out[k += 1] = (x[i] + x[j]) / 2
    end
    return sort!(out)
end

"""
    cross_differences(x, y) -> Vector

The ``n_x n_y`` differences ``xᵢ - yⱼ``, sorted ascending. Their median is the
two-sample Hodges-Lehmann estimator, and the distribution-free interval for the
shift is a pair of their order statistics.
"""
function cross_differences(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real,S<:Real}
    nx, ny = length(x), length(y)
    m = nx * ny
    check_estimate_count(m, "cross-group differences")
    out = Vector{float(promote_type(T, S))}(undef, m)
    k = 0
    @inbounds for xi in x, yj in y
        out[k += 1] = xi - yj
    end
    return sort!(out)
end

## INDEX RULES
#
# Both intervals are read off a sorted set of pairwise estimates `vals` of length `m` as
# `(vals[k + 1], vals[m - k])`. All that differs between the exact and the
# approximate test is how `k` is chosen.

# Largest `k` in `0:H` at which `pred` still holds, for `pred` true on an initial
# segment. Binary search, because each `pred` call runs a whole lattice DP.
function last_true(H::Integer, pred)
    pred(0) || return 0
    lo, hi = 0, Int(H)
    while lo < hi
        mid = lo + div(hi - lo + 1, 2)
        if pred(mid)
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

"""
    exact_ci_index(m, alpha, cdf) -> Int

Invert the exact null distribution conservatively: the largest `k` whose interval
still attains at least `1 - alpha`. `cdf(k)` is the null probability of a statistic
at most `k`.

This is the rule R's `wilcox.test` uses, where the lower endpoint is
`diffs[qsignrank(alpha/2, n)]`; taking `k = qsignrank(alpha/2, n) - 1` reproduces it.
"""
function exact_ci_index(m::Integer, alpha::Real, cdf)
    half = alpha / 2
    return last_true(div(m, 2), k -> cdf(k) < half)
end

"""
    normal_ci_index(m, mu, sigma, alpha) -> Int

Invert the normal approximation to the null distribution, `mu` and `sigma` being its
mean and (tie-corrected) standard deviation.

The exact rule takes the lower endpoint to be the `Cα`-th order statistic, for
`Cα = min{j : P(W ≤ j) ≥ alpha/2}`. The statistic is supported on a lattice, so
`P(W ≤ j) ≈ Φ((j + 1/2 - mu) / sigma)` and that critical value approximates to

    Cα = ceil(mu - z * sigma - 1/2)

with `k = Cα - 1`. This is the continuity correction R's `wilcox.test` applies by
default. Dropping it makes the interval narrower than the exact one for 7 of 66 signed
rank sample sizes at `level = 0.95` and 45 of 66 at `level = 0.90`; with it the interval
is never narrower than exact at `level = 0.95` or above.
"""
function normal_ci_index(m::Integer, mu::Real, sigma::Real, alpha::Real)
    z = quantile(Normal(), 1 - alpha / 2)
    return clamp(ceil(Int, mu - z * sigma - 0.5) - 1, 0, div(m, 2))
end

# `level`/`tail` to the two-sided alpha the index rules work in. A one-sided bound
# at `level` is the corresponding endpoint of the two-sided interval at `2 level - 1`.
function ci_alpha(level::Real, tail::Symbol)
    check_level(level)
    check_tail(tail)
    return tail === :both ? 1 - level : 2 * (1 - level)
end

function ci_from_estimates(vals::AbstractVector, k::Integer, tail::Symbol)
    m = length(vals)
    left = vals[k + 1]
    right = vals[m - k]
    if tail === :both
        return (left, right)
    elseif tail === :left
        return (left, oftype(left, Inf))
    else # tail === :right
        return (oftype(right, -Inf), right)
    end
end

## HODGES-LEHMANN

"""
    hodgeslehmann(test)

The Hodges-Lehmann estimator underlying `test`: the median of the Walsh averages for
the signed rank tests, and the median of the cross-group differences for the
Mann-Whitney U tests.

This is the point estimate the distribution-free interval returned by
[`confint`](@ref) is built around, and the `estimate` R's `wilcox.test` reports where it
takes its exact route. On its approximate route R solves for that estimate numerically
instead, and reports a slightly different number: `9.71184` against `9.675` on the
fifteen-point sample in the tests.

It is generally not the sample median.
"""
function hodgeslehmann end
