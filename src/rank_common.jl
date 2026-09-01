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
    elseif applicable(method, stats)
        # `applicable` rather than `method isa Function`, which asks about type ancestry
        # and not about callability: a struct with a call method is callable but is not
        # under `Function` unless it says so, and a parameterised rule is the reason this
        # keyword exists. The symbol branches are tested first, so nothing that should be
        # a symbol reaches here, and a callable whose signature does not match falls
        # through to the message below rather than raising a MethodError from in here.
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
unbounded time, and `method = :approximate` is the way on.

`confint` does not enumerate, since it inverts the tie-free distribution whichever route the
p-value took, so this bound does not reach it. Its own cost is bounded by
[`MAX_EXACT_CI_ESTIMATES`](@ref). See #7.

At the bound itself a tied signed rank p-value takes about 14 s over its ``2^{30}`` sign
assignments, and the largest balanced split it admits, 16 against 16, about 20 s over
`binomial(32, 16)` rank assignments. One observation further doubles both.
"""
const MAX_EXACT_ENUMERATION_N = 30

function check_exact_enumeration(n::Integer)
    n <= MAX_EXACT_ENUMERATION_N && return nothing
    throw(ComputationTooLarge(
        "an exact test with ties enumerates all 2^n sign assignments, which is not " *
        "feasible for n = $n non-zero observations (limit " *
        "MAX_EXACT_ENUMERATION_N = $MAX_EXACT_ENUMERATION_N). Pass `method = :approximate`."))
end

function check_exact_enumeration(nx::Integer, ny::Integer)
    count = rank_binomial(nx, ny)
    (count !== nothing && count <= Int128(2)^MAX_EXACT_ENUMERATION_N) && return nothing
    throw(ComputationTooLarge(
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

The bound is on memory alone, which is what the approximate routes spend: forming this many
estimates costs about 8 MB, and sorting them about that again. It admits a signed rank
sample of roughly 1400 observations and a two-sample design of 1000 against 1000. The
two-sample exact interval spends time as well and carries its own bound below this,
[`MAX_EXACT_CI_ESTIMATES`](@ref); the one-sample exact interval bisects cheaply enough
that this bound is the only one it meets, about 6 s at the boundary.
"""
const MAX_PAIRWISE_ESTIMATES = 1_000_000

function check_estimate_count(m::Integer, what::AbstractString)
    m <= MAX_PAIRWISE_ESTIMATES && return nothing
    throw(ComputationTooLarge(
        "refusing to form $m $what: more than MAX_PAIRWISE_ESTIMATES = $MAX_PAIRWISE_ESTIMATES. " *
        "The distribution-free interval and the Hodges-Lehmann estimator are computed " *
        "by materialising this set, which is quadratic in the sample size (see #7)."))
end

"""
Largest set of pairwise estimates whose *exact* Mann-Whitney interval will be inverted.

Selecting an endpoint from the exact null distribution runs a lattice recursion per
bisection step, and the two-sample recursion (`wilcoxcdf`) is the expensive one: at this
bound, a 300-against-300 design, the interval takes about 11 s, from 0.6 s at 158 against
158, and it grows steeply past it, 35 s at 400 against 400. Past it the exact interval is
refused; the approximate interval, which inverts a closed form, is not.

The one-sample recursion (`signrankcdf`) is cheap enough that the signed rank interval
carries no time bound at all: it reaches the memory bound
[`MAX_PAIRWISE_ESTIMATES`](@ref), a sample of 1413, in about 6 s.

The number was first 25\\_000, sized against the ``m/2``-recursion endpoint scan that
inverted the signed rank interval before #361; the bisection that replaced the scan
moved the cost from the one-sample route to the two-sample one, and the same wall-clock
ceiling, re-measured, sits at 90\\_000.
"""
const MAX_EXACT_CI_ESTIMATES = 90_000

function check_exact_ci_cost(m::Integer, escape::AbstractString)
    m <= MAX_EXACT_CI_ESTIMATES && return nothing
    throw(ComputationTooLarge(
        "refusing to invert the exact null distribution over $m pairwise estimates: more " *
        "than MAX_EXACT_CI_ESTIMATES = $MAX_EXACT_CI_ESTIMATES. Choosing an endpoint runs " *
        "a two-sample lattice recursion per bisection step, so this route is bounded below " *
        "MAX_PAIRWISE_ESTIMATES = $MAX_PAIRWISE_ESTIMATES, which bounds memory alone. " *
        escape))
end

"""
    walsh_averages(x) -> Vector

The ``n(n+1)/2`` Walsh averages ``(xᵢ + xⱼ)/2`` for ``i ≤ j``, sorted ascending.
Their median is the one-sample Hodges-Lehmann estimator, and the distribution-free
interval for the pseudomedian is a pair of their order statistics.
"""
function walsh_averages(x::AbstractVector{T}) where T<:Real
    # the triangular loop below indexes `x` over `1:n` with bounds checking off, so state
    # the assumption rather than corrupting memory if it is ever broken. `cross_differences`
    # needs no such guard: it iterates the values and indexes only its own fresh output.
    Base.require_one_based_indexing(x)
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
    exact_ci_index(m, alpha, cdf; tail = :both) -> Int

Invert the exact null distribution conservatively: the largest `k` whose interval
still attains at least `1 - alpha`. `cdf(k)` is the null probability of a statistic
at most `k`.

This is the rule R's `wilcox.test` uses, where the lower endpoint is
`diffs[qsignrank(alpha/2, n)]`; taking `k = qsignrank(alpha/2, n) - 1` reproduces it.

Small samples cannot reach every level: the widest interval this form admits,
`k = 0`, has coverage `1 - 2 cdf(0)`, which is `1 - 2^(1-n)` for an untied signed rank
sample, so `0.95` is out of reach below `n = 6` and `0.99` below `n = 8`. That interval
is still the best available and is returned, with a warning naming the coverage it
actually has, rather than being returned as though it met the request. R warns here too.

`tail` only shapes that warning: `alpha` is already the two-sided mass whichever tail
the interval keeps, but the level a one-sided caller asked for is `1 - alpha/2`, not
`1 - alpha`, and the coverage its single overshooting tail costs is `cdf(0)` once,
not twice.
"""
function exact_ci_index(m::Integer, alpha::Real, cdf; tail::Symbol = :both)
    half = alpha / 2
    p0 = cdf(0)
    if p0 > half
        # `1 - 2 p0` double-counts the atom at 0 when the null distribution is degenerate,
        # so it can come out negative; the coverage is nowhere below zero.
        sides = tail === :both ? 2 : 1
        @warn "the requested confidence level is not attainable for a sample this small; " *
              "returning the widest interval this construction admits" requested = 1 - sides * half attainable =
              max(0.0, 1 - sides * p0)
        return 0
    end
    return last_true(div(m - 1, 2), k -> cdf(k) < half)
end

"""
    normal_ci_index(m, mu, sigma, alpha; tail = :both) -> Int

Invert the normal approximation to the null distribution, `mu` and `sigma` being its
mean and (tie-corrected) standard deviation.

The exact rule takes the lower endpoint to be the `Cα`-th order statistic, for
`Cα = min{j : P(W ≤ j) ≥ alpha/2}`. The statistic is supported on a lattice, so
`P(W ≤ j) ≈ Φ((j + 1/2 - mu) / sigma)` and that critical value approximates to

    Cα = ceil(mu - z * sigma - 1/2)

with `k = Cα - 1`. `mu - z * sigma` is the `alpha/2` quantile of that normal, which is
how it is read below, off the parameterised distribution rather than by scaling the
standard one. This is the continuity correction R's `wilcox.test` applies by
default. Dropping it makes the interval narrower than the exact one for 7 of 66 signed
rank sample sizes at `level = 0.95` and 45 of 66 at `level = 0.90`; with it the interval
is never narrower than exact at `level = 0.95` or above.

Below that it can still be narrower by one order statistic, on strongly unbalanced
designs: at `level = 0.90` that happens for 131 of the 1444 Mann-Whitney shapes with
`nx, ny` in `3:40`, concentrated at `nx = 3`. It is the exact rule that is conservative
there rather than this one that is loose, and coverage stays at nominal rather than
falling below it: at `(3, 12)`, `(3, 15)` and `(3, 21)` the normal route covers `0.902`,
`0.901` and `0.900` against the exact route's `0.932`, `0.927` and `0.919`.

A negative `k` means the approximation puts the endpoint outside the pairwise estimates, so
the widest interval available is returned instead, with a warning. The warning does not
claim the level is unattainable, because that does not follow: at `n = 8` and `level = 0.99` this
rule asks for `k = -1` while the exact route reaches `0.9922`, which meets the request.
What it does mean is that the sample is too small for this route, and `method = :exact` is
where to go.

As in [`exact_ci_index`](@ref), `tail` only shapes the warning: a one-sided caller asked
for level `1 - alpha/2`, and that is what the warning should name as the request.
"""
function normal_ci_index(m::Integer, mu::Real, sigma::Real, alpha::Real; tail::Symbol = :both)
    q = quantile(Normal(mu, sigma), alpha / 2)
    k = ceil(Int, q - one(q) / 2) - 1
    if k < 0
        sides = tail === :both ? 2 : 1
        @warn "the normal approximation puts the interval endpoint outside the set of " *
              "pairwise estimates for a sample this small; returning the widest interval " *
              "it admits, which may not reach the requested level. The exact test is " *
              "the way on" requested = 1 - sides * (alpha / 2)
    end
    return clamp(k, 0, div(m - 1, 2))
end

# `level`/`tail` to the two-sided alpha the index rules work in. A one-sided bound
# at `level` is the corresponding endpoint of the two-sided interval at `2 level - 1`.
function ci_alpha(level::Real, tail::Symbol)
    # the rank tests are the only `confint` methods here declaring `level::Real`; every
    # other one declares `Float64` and gets this conversion for free at its own boundary,
    # so `check_level` stays `Float64`-only and the conversion happens here instead
    lev = Float64(level)
    check_level(lev)
    check_tail(tail)
    return tail === :both ? 1 - lev : 2 * (1 - lev)
end

# A one-sided interval keeps the endpoint that inverts the test of the same name, which is
# the opposite endpoint to the one the tail is named after: the alternative `tail = :left`
# stands for is location *below* the null, and what is compatible with that is an upper
# bound. Eight of the nine other `confint` methods here that take a `tail` do the same, as
# does R's `wilcox.test` under `alternative = "less"`. These four returned the other
# endpoint until #368. `SignTest` is the ninth and still returns a lower bound, capped at
# the null rather than infinite; that is the same defect and is filed separately.
function ci_from_estimates(vals::AbstractVector, k::Integer, tail::Symbol)
    m = length(vals)
    left = vals[k + 1]
    right = vals[m - k]
    if tail === :both
        return (left, right)
    elseif tail === :left
        return (oftype(right, -Inf), right)
    else # tail === :right
        return (left, oftype(left, Inf))
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
