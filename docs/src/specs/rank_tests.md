# Rank-based location inference

Specification of two rank procedures and the location inference built on them:

  - the **one-sample procedure**, the Wilcoxon signed rank test, for a single sample of
    differences, and equally for **paired** data, which enters as the within-pair
    differences ``d_i = x_i - y_i``;
  - the **two-sample procedure**, the Wilcoxon rank sum test, equivalently the
    Mann-Whitney ``U`` test, for two independent samples.

[§2](@ref "2. The one-sample procedure (Wilcoxon signed rank)") and
[§3](@ref "3. The two-sample procedure (Wilcoxon rank sum, Mann-Whitney U)") take the two
procedures in turn, each giving the statistic, its exact null distribution together with
the two forms in which it is computed and the normal approximation to it, and the p-value
read off whichever is in force. The location inference is then common to both: the
pairwise estimates it is read off ([§4](@ref "4. Pairwise estimates")), the Hodges-Lehmann point
estimate ([§5](@ref "5. Point estimation")), and the confidence interval obtained by
inverting either test ([§6](@ref "6. Interval estimation")). The rest relates the
specification to other implementations
([§7](@ref "7. Relation to other implementations")) and gives worked values for checking
one ([§8](@ref "8. Worked values for the rank tests")).

**In this package.** The one-sample procedure is [`ExactSignedRankTest`](@ref) and
[`ApproximateSignedRankTest`](@ref), the two-sample procedure
[`ExactMannWhitneyUTest`](@ref) and [`ApproximateMannWhitneyUTest`](@ref);
[`SignedRankTest`](@ref) and [`MannWhitneyUTest`](@ref) are functions that pick one of those
four rather than types of their own. [§9](@ref "9. The rank tests in this package") maps the
specification onto them and records where they depart from it.

## 1. Preliminaries and notation

| symbol | meaning |
|:---|:---|
| ``D_1, \dots, D_N`` | the one-sample observations as random variables |
| ``d_1, \dots, d_N`` | the values they took, which is what is computed from; for paired data ``d_i = x_i - y_i`` |
| ``n`` | the number of ``d_i`` with ``d_i \ne 0`` |
| ``X_1, \dots, X_{n_x}``, ``Y_1, \dots, Y_{n_y}`` | the two-sample observations as random variables |
| ``x_1, \dots, x_{n_x}``, ``y_1, \dots, y_{n_y}`` | the values they took |
| ``N`` | the total number of observations: the number of ``d_i``, zeros included, one-sample; ``n_x + n_y`` two-sample |
| ``\alpha`` | two-sided error rate; the interval has nominal coverage ``1-\alpha`` |
| ``z_q`` | the ``q`` quantile of the standard normal distribution |
| ``\Phi`` | the standard normal CDF |
| ``\ell`` | the number of values being ranked at once |
| ``m`` | the number of pairwise estimates of [§4](@ref "4. Pairwise estimates") |
| ``V_{(1)} \le \dots \le V_{(m)}`` | those pairwise estimates, sorted |
| ``A_{ij}``, ``D_{ij}`` | the pairwise estimates before sorting ([§4.1](@ref "4.1 Definitions")); the two indices distinguish ``D_{ij}`` from the random variable ``D_i`` |

**Random variables and their values.** Capitals denote random variables and lowercase the
values they took: ``D_i`` is an observation of the one-sample procedure regarded as random,
``d_i`` the number in hand. Every hypothesis and every null distribution on this page is a
statement about the capitals; everything computed, the ranks, the statistics, the estimates
and the intervals, is a function of the lowercase. The statistics ``W^+`` and ``U`` are
written as capitals throughout, since which is meant is fixed by context: inside a
probability they are random, and elsewhere they are the observed value.

**Ranks and midranks.** Ranking always happens on a single vector of ``\ell`` values:
the ``\ell = n`` retained ``|d_i|`` in
[§2](@ref "2. The one-sample procedure (Wilcoxon signed rank)"), and the ``\ell = N``
pooled observations in
[§3](@ref "3. The two-sample procedure (Wilcoxon rank sum, Mann-Whitney U)").

Sort those values into ascending order. Absent ties, the **rank** of a value is its
position in that order: the sort gives each of the ``\ell`` values a distinct position in
``1, \dots, \ell``, so ranking is a bijection from the values onto those positions and the
ranks it produces are a permutation of them. Tied values occupy a block of
consecutive positions, and each is instead given the average of the positions in its block,
its **midrank**. For a vector ``v``,

```math
R_i = \#\{ j : v_j < v_i \} + \frac{1 + \#\{ j : v_j = v_i \}}{2} ,
```

the second count including ``i`` itself. A value tied with no other therefore has
``\#\{ j : v_j = v_i \} = 1`` and receives its ordinary rank, while a group of ``t``
equal values occupying positions ``r+1, \dots, r+t`` receives the common midrank
``r + (t+1)/2``. The two notions differ only on tied data.

**Tie totals.** For a vector ``v``, write

```math
T(v) = \sum_{g} (t_g^3 - t_g)
```

where ``t_g`` is the multiplicity of the ``g``-th group of equal values in ``v``. Since
``t^3 - t = t(t-1)(t+1)``, groups of size one contribute nothing, and ``T(v) = 0`` exactly
when ``v`` has no ties. Both normal approximations reduce the tie pattern to
``T`` and use nothing else from it, reaching it through the null variance of their
statistic ([§2.1](@ref "2.1 Model, estimand, statistic"),
[§3.1](@ref "3.1 Model, estimand, statistic")). Neither exact distribution uses ``T``: under
ties the permutation one conditions on the whole multiset of midranks and enumerates
([§2.2.2](@ref "2.2.2 The permutation distribution"),
[§3.2.2](@ref "3.2.2 The permutation distribution")).

### 1.1 Mathematical observations

Three properties of midranks, used repeatedly below.

**They depend only on the multiset of values.** The closed form counts values and nothing
else, so no input ordering and no tie-breaking convention enters: permuting ``v`` permutes
``R`` the same way, and equal values always receive equal ranks.

**Ties leave the rank total alone.** Whether or not ``v`` has ties,
``\sum_i R_i = \ell(\ell+1)/2``, since a group of tied values contributes one copy of the
average of the positions it occupies for each position it occupies. This is why ties never
affect the null mean of either statistic.

**Ties shrink the rank spread by exactly ``T(v)/12``.** Consider a group of ``t`` tied
values occupying positions ``r+1, \dots, r+t``. Untied, they would have carried those ``t``
distinct positions as their ranks, contributing ``\sum_{k=1}^{t} (r+k)^2`` to the sum of
squares ``\sum_i R_i^2``; tied, each carries the common midrank ``r + (t+1)/2`` instead,
contributing ``t`` copies of its square. Replacing distinct values by their average always
lowers a sum of squares, here by

```math
\sum_{k=1}^{t} (r+k)^2 \;-\; t\left(r + \tfrac{t+1}{2}\right)^{2}
  = \frac{t(t+1)(2t+1)}{6} - \frac{t(t+1)^2}{4}
  = \frac{t^3 - t}{12} ,
```

independent of ``r``: only the size of a group matters, not where it sits. Summing over
groups,

```math
\sum_i R_i^2 = \frac{\ell(\ell+1)(2\ell+1)}{6} - \frac{T(v)}{12} ,
\qquad
\sum_i (R_i - \bar R)^2 = \frac{\ell(\ell^2-1)}{12} - \frac{T(v)}{12} .
```

So the cube in ``T`` is not a fitted constant, and ``T`` is not a correction bolted onto
the variances: it is twelve times the amount by which ties shrink the spread of the ranks.
Every appearance of ``T`` below is one of these two identities substituted into a variance,
which is why it always arrives divided by ``12``, or by ``48 = 4 \times 12``.

### 1.2 A worked ranking

Both procedures reduce their data to midranks before anything else happens, and both
summarise the ties through ``T``. A small example of each fixes the mechanics. The two
statistics are defined here as they are computed;
[§2.1](@ref "2.1 Model, estimand, statistic") and
[§3.1](@ref "3.1 Model, estimand, statistic") restate them alongside their supports, null
moments and null distributions.

**One sample.** Take ``d = (2.1,\, -0.7,\, 0,\, 1.4,\, -2.1,\, 0.7)``. The zero is
discarded, leaving ``n = 5``, and the remaining absolute values are ranked. Two pairs tie:
the two ``0.7`` occupy positions ``1`` and ``2``, the two ``2.1`` occupy ``4`` and ``5``,
and each takes the average of the positions it spans.

| ``d_i`` | ``2.1`` | ``-0.7`` | ``0`` | ``1.4`` | ``-2.1`` | ``0.7`` |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| ``\lvert d_i \rvert`` | ``2.1`` | ``0.7`` | discarded | ``1.4`` | ``2.1`` | ``0.7`` |
| midrank ``R_i`` | ``4.5`` | ``1.5`` | | ``3`` | ``4.5`` | ``1.5`` |
| counted in ``W^+`` | yes | no | | yes | no | yes |

The two tied pairs give ``T(\lvert d \rvert) = (2^3 - 2) + (2^3 - 2) = 12``.

The **signed rank statistic** adds up the midranks carried by the positive observations,
discarding the rest:

```math
W^+ = \sum_{i \,:\, d_i > 0} R_i .
```

Here that is the three cells marked *yes*, carrying midranks ``4.5``, ``3`` and ``1.5``,
so ``W^+ = 9``. Its null mean and variance
([§2.1](@ref "2.1 Model, estimand, statistic")) come to ``7.5`` and ``13.5``, the variance
having been reduced from its untied value ``13.75`` by ``T(|d|)/48``.

**Two samples.** Take ``x = (4.2,\, 1.5,\, 3.7,\, 5.4,\, 2.2)`` and
``y = (6.1,\, 2.8,\, 3.7,\, 4.9)``, so ``n_x = 5``, ``n_y = 4`` and ``N = 9``. Here the
two samples are pooled before ranking, and the tie falls *across* them: the two ``3.7``
share positions ``4`` and ``5``.

| pooled value | ``1.5`` | ``2.2`` | ``2.8`` | ``3.7`` | ``3.7`` | ``4.2`` | ``4.9`` | ``5.4`` | ``6.1`` |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| from | ``x`` | ``x`` | ``y`` | ``x`` | ``y`` | ``x`` | ``y`` | ``x`` | ``y`` |
| midrank | ``1`` | ``2`` | ``3`` | ``4.5`` | ``4.5`` | ``6`` | ``7`` | ``8`` | ``9`` |

The **Mann-Whitney statistic** adds up the pooled midranks that fell to ``x``, then
subtracts the smallest total those ``n_x`` ranks could possibly have taken:

```math
U = \sum_{i=1}^{n_x} R_i - \frac{n_x(n_x+1)}{2} ,
```

with ``R_i`` the pooled midrank of ``x_i``. Here ``x`` received
``\{1,\, 2,\, 4.5,\, 6,\, 8\}``, totalling ``21.5``, and the subtracted minimum is
``1 + 2 + 3 + 4 + 5 = 15``, fixed by ``n_x`` alone, so ``U = 6.5``. That subtraction turns
the rank sum into a count, and counting directly over the ``n_x n_y = 20`` pairs gives the
same number: ``5.4`` beats ``2.8``, ``3.7`` and ``4.9``; ``4.2`` beats ``2.8`` and ``3.7``;
``3.7`` beats ``2.8`` and ties the other ``3.7``; for ``6 + \tfrac{1}{2} = 6.5``. That the
rank sum and the pair count agree is not a coincidence of this example:
[§3.1](@ref "3.1 Model, estimand, statistic") proves they are equal for every sample.
One tied group of two gives ``T([x; y]) = 6``, and the null mean and variance are ``10``
and ``16.52778``.

Both examples are tied, which is the case worth working. Absent ties every group has size
one, ``T = 0``, the midranks are the ordinary ranks ``1, \dots, \ell``, and both variances
collapse to their first term: ``13.75`` here in place of ``13.5``, and ``16.66667`` in
place of ``16.52778``. The correction is small on samples this size and stays small, but it
is the only route by which the tie pattern reaches either normal approximation.

## 2. The one-sample procedure (Wilcoxon signed rank)

The signed rank statistic uses both the sign of each ``d_i`` and the rank of its magnitude
among the others, which is what a single sample has to offer.
[§2.1](@ref "2.1 Model, estimand, statistic") fixes the model, the estimand and the
statistic, [§2.2](@ref "2.2 Null distribution of the signed rank statistic") gives the
distribution that statistic has under the null, and [§2.3](@ref "2.3 p-values") reads the
p-values off it.
[§3](@ref "3. The two-sample procedure (Wilcoxon rank sum, Mann-Whitney U)") follows the
same order for the two-sample statistic.

### 2.1 Model, estimand, statistic

**Model.** The observations are independent random variables ``D_1, \dots, D_N``. The null
hypothesis is that the distribution of each is symmetric about ``0``, meaning
``\mathbb{P}(D_i > t) = \mathbb{P}(D_i < -t)`` for every ``t``. What that requires, and what it does not, is
[§2.1.1](@ref "2.1.1 Assumptions").

**Estimand.** The **pseudomedian** of ``F``: the median of ``(D + D')/2`` for independent
``D, D' \sim F``. If ``F`` is symmetric about ``\theta`` then the pseudomedian is
``\theta``, and coincides with the median of ``F`` and, where it exists, with the mean.
For asymmetric ``F`` the pseudomedian is a different functional from the median, and it is
the pseudomedian that [§5](@ref "5. Point estimation") and
[§6](@ref "6. Interval estimation") estimate.

**Hypotheses.** The null the level is exact for is the one in the Model paragraph above, a
statement about the distributions of the ``D_i`` rather than about a parameter of them:

```math
H_0 : \mathbb{P}(D_i > t) = \mathbb{P}(D_i < -t) \quad \text{for every } t \text{ and every } i .
```

Under the location model ``F(t) = F_0(t - \theta)`` with ``F_0`` symmetric about ``0``, so
that ``\theta`` is the pseudomedian, that null reads ``\theta = 0`` and the three
alternatives are

```math
H_1 : \theta \ne 0 \quad (\texttt{tail = :both}), \qquad
H_1 : \theta < 0 \quad (\texttt{tail = :left}), \qquad
H_1 : \theta > 0 \quad (\texttt{tail = :right}) ,
```

small ``W^+`` being the evidence for ``\theta < 0``. The same names carry over to the
one-sided intervals of [§6.5](@ref "6.5 One-sided intervals"), which keep the endpoint that
inverts the test of the same name: ``\theta < 0`` admits an upper bound, so `tail = :left`
returns one.

To test ``\theta = \theta_0`` for a value other than zero, run the procedure on
``d_i - \theta_0``. Nothing here provides for it directly, and nothing needs to: every
quantity below is computed from the sample it is handed, and
[§4.2](@ref "4.2 The counting identity") is that shift written out.

**Zeros.** Observations with ``d_i = 0`` are discarded before ranking; a zero has no sign
to rank. All of [§2](@ref "2. The one-sample procedure (Wilcoxon signed rank)") and all of
[§4](@ref "4. Pairwise estimates")–[§6](@ref "6. Interval estimation") are computed from the
``n`` retained observations. ``N`` enters nothing this specification computes, so an
implementation reporting a sample size should report both counts, or say which one it means.
(This package reports both, and does compute one thing from all ``N``: see
[§9](@ref "9. The rank tests in this package").)

The discard is total. The procedure returns exactly what it would have returned had those
observations never been recorded, so appending zeros to a sample moves neither ``W^+``, nor
the p-value, nor the estimate, nor the interval.

Discarding conditions the test on which observations were non-zero, and so on ``n``, in
the sense set out in [§2.2.2](@ref "2.2.2 The permutation distribution"). That costs
nothing in exactness, since the retained signs are still independent and each ``\pm`` with
probability ``1/2`` whatever the zeros did. It does cost information: only ``n`` of the
``N`` observations reach the test, and the support of a statistic built from ``n`` ranks is
coarser than one built from ``N``, so the smallest attainable p-value is larger and the
power lower.

That is a choice rather than a consequence, and it is a strong one, since a difference of
exactly ``0`` is the observation most consistent with symmetry about ``0`` and is given no
weight at all. The alternative is Pratt's convention [pratt1959](@cite), which ranks the
zeros along with the rest and then omits their ranks from the sum: the zeros inflate the
midranks of the retained observations and so do reach the p-value. The two conventions
disagree often on samples containing zeros, and in both directions, so neither is uniformly
the more conservative. This page specifies the discard, which is Wilcoxon's own convention
and the one R's `wilcox.test` follows.

**Statistic.** Re-index the retained observations as ``d_1, \dots, d_n`` and let ``R_i``
be the midrank of ``|d_i|`` among those ``n`` absolute values. Then

```math
W^+ = \sum_{\substack{1 \le i \le n \\ d_i > 0}} R_i ,
```

the index running over the retained observations only. Both the ranking and the sum are
taken after the zeros are gone, so a zero is not merely left out of the sum: it is absent
from the ranking the sum is over, which is what makes the discard above total. Under
Pratt's convention it would be the other way round, present in the ranking and left out of
the sum.

It runs from ``0``, every retained observation negative, to ``n(n+1)/2``, every one
positive, and is symmetric about ``n(n+1)/4`` under the null. Absent ties it takes integer
values, so the support is ``\{0, 1, \dots, n(n+1)/2\}``; under ties the midranks are
half-integers and so is ``W^+``. At ``n = 2`` with the two absolute values equal, both
midranks are ``1.5`` and the support is ``\{0,\, 1.5,\, 3\}``.

Under the null the signs are independent of ``|D|`` and each ``\pm`` with probability
``1/2``, so conditionally on the midranks ``W^+ = \sum_i R_i B_i`` with
``B_i \overset{\text{iid}}{\sim} \mathrm{Bernoulli}(1/2)``. Hence
``\mathbb{E}[W^+] = \tfrac{1}{2}\sum_i R_i`` and
``\operatorname{Var}(W^+) = \tfrac{1}{4}\sum_i R_i^2``, and substituting
[§1.1](@ref "1.1 Mathematical observations"),

```math
\mathbb{E}[W^+] = \frac{n(n+1)}{4}, \qquad
\operatorname{Var}(W^+) = \frac{n(n+1)(2n+1)}{24} - \frac{T(|d|)}{48} .
```

The mean is untouched by ties because the rank total is; the variance is not, because the
sum of squares is not.

#### 2.1.1 Assumptions

Symmetry about ``0`` is the whole of what
[§2.2](@ref "2.2 Null distribution of the signed rank statistic") and
[§2.3](@ref "2.3 p-values") use: it makes the signs independent of ``|D|`` and each ``\pm``
with probability ``1/2``, and every null distribution below follows from that alone. In
particular the ``D_i`` need not be identically distributed for the null distribution to
hold, though [§5](@ref "5. Point estimation") and [§6](@ref "6. Interval estimation") do
assume a common ``F``, since the pseudomedian they estimate is a functional of one
distribution.

Continuity of ``F`` is a convenience rather than a requirement. It makes zeros and ties
events of probability zero, so the ranks are the integers ``1, \dots, n`` and the lattice
distribution of [§2.2.1](@ref "2.2.1 The lattice distribution") applies. Under a
discrete or mixed ``F`` the procedure remains exact provided ties are handled as in
[§2.2.2](@ref "2.2.2 The permutation distribution"), which also explains in what sense
such a test is exact, and zeros are discarded as in the **Zeros** paragraph above.

What cannot be weakened is symmetry. Testing ``\operatorname{median}(F) = 0`` is a
different problem, and this test does not solve it. ``W^+`` sits at ``n(n+1)/4`` on average
only when ``F`` is symmetric; an asymmetric ``F`` with median ``0`` has a pseudomedian away
from ``0``, and ``W^+`` follows the pseudomedian, so it is centred somewhere else and the
test rejects more often than ``\alpha``. Take ``F`` the law of
``\mathrm{Exponential}(1) - \ln 2``, whose median is ``0`` and whose pseudomedian is about
``0.146``: at ``n = 20`` the nominal ``0.05`` two-sided test rejects about ``10\%`` of the
time.

### 2.2 Null distribution of the signed rank statistic

A p-value needs the distribution of ``W^+`` under ``H_0``. For a given sample there is
exactly one exact null distribution; what varies is the form in which it can be computed,
and whether it is computed at all rather than approximated. The three forms carry names,
used in place of section numbers from this point on:

| name | what it is | available for | cost |
|:---|:---|:---|:---|
| **lattice** ([§2.2.1](@ref "2.2.1 The lattice distribution")) | the exact distribution in its untied form, tabulated by a recursion | untied samples only, ``T(\lvert d \rvert) = 0`` | ``O(n^3)`` |
| **permutation** ([§2.2.2](@ref "2.2.2 The permutation distribution")) | the exact distribution in general, enumerated from the sample | any sample, and the only exact form once ``T(\lvert d \rvert) > 0`` | ``O(2^n)`` |
| **normal** ([§2.2.3](@ref "2.2.3 The normal approximation")) | an approximation to it, from the moments above | any sample | ``O(n)`` |

The tie pattern decides which forms are *available*; it does not by itself decide which is
used. The costs are per p-value, once the sample is ranked.

Lattice and permutation are two computations of one distribution, not two distributions,
which is why the table says form and not choice of null. The permutation form conditions on
the multiset of midranks the sample produced; absent ties that multiset is
``\{1, \dots, n\}`` for *every* sample of size ``n``, so conditioning on it fixes
nothing, the distribution is unconditional, and it depends on the sample through ``n``
alone. That special case is the lattice form, tabulated once where the permutation form is
rebuilt for each sample. On untied data the enumeration of
[§2.2.2](@ref "2.2.2 The permutation distribution") reproduces the recursion of
[§2.2.1](@ref "2.2.1 The lattice distribution") to the last digit, at far greater cost, which
is why the recursion is taken whenever it applies.

**Choosing between them.** There are two decisions here, and they are of different kinds.

**Lattice against permutation** is forced by the data and is not a choice of null: with no
ties the midranks are ``1, \dots, n`` and the recursion applies, and with ties they are
not, so it does not. The two forms agree wherever both apply, and
[§2.3](@ref "2.3 p-values") reads the same three formulas off either.

**Exact against normal** *is* a choice, and it is where the sample size enters: it trades
the cost above, prohibitive for the enumeration on a large tied sample, against a p-value
that is only asymptotically right. An untied sample of any size may still be taken through
the normal, and a large one usually is. The normal is also the only form that
responds to ties through ``T``, since the exact distribution conditions on the midranks
rather than summarising them.

**In this package.** Two types implement the procedure, `ExactSignedRankTest` and
`ApproximateSignedRankTest`. `SignedRankTest` is not a third: it is a function that ranks
the sample, applies the rule below, and returns one of those two. The tie pattern and the
number of non-zero differences are the only properties of the data that enter:

| call | ``\lvert d \rvert`` untied | ``\lvert d \rvert`` tied |
|:---|:---|:---|
| `SignedRankTest(d)` | lattice for ``n \le 50``, normal above | permutation for ``n \le 15``, normal above |
| `ExactSignedRankTest(d)` | lattice, no bound | permutation, up to ``n = 30`` |
| `ApproximateSignedRankTest(d)` | normal | normal |

Reading the first row: the normal approximation is what `SignedRankTest` falls back on when
the sample is too large for the exact route it would otherwise take, and ties lower that size
limit from ``50`` to ``15``, because they replace the lattice recursion with an enumeration
costing ``2^n``. A sample of ``40`` differences therefore goes exact if untied and normal if
any two of the ``\lvert d_i \rvert`` are equal. Here ``n`` counts the non-zero differences, so
the zeros of [§2.1](@ref "2.1 Model, estimand, statistic") are discarded before the rule is
applied.

Because the rule reads the data, a `SignedRankTest` call has no one return type: the same
expression yields an `ExactSignedRankTest` on one sample and an `ApproximateSignedRankTest`
on another. `method = :exact` or `method = :approximate` overrides the rule and fixes which,
which is what an analysis that must reproduce across versions should do; naming the type
directly does the same. A callable, passed `(; n, n_nonzero, ties, tie_adjustment)` and
returning one of those two symbols, replaces the rule with one of the caller's own.

Forcing the exact route on a large tied sample is declined rather than served slowly. Beyond
[`MAX_EXACT_ENUMERATION_N`](@ref HypothesisTests.MAX_EXACT_ENUMERATION_N) ``= 30`` non-zero
differences, `pvalue` on a tied `ExactSignedRankTest` throws
[`ComputationTooLarge`](@ref HypothesisTests.ComputationTooLarge) naming ``n``
and pointing at `method = :approximate`, rather than starting an enumeration whose cost
doubles with every further observation: at ``n = 40`` there would be
``2^{40} \approx 1.1 \times 10^{12}`` sign patterns, and the call would not return. The
bound is a property of this enumeration rather than of the problem;
[§2.2.2](@ref "2.2.2 The permutation distribution") closes with the polynomial algorithm
that would remove it, proposed as
[issue #370](https://github.com/JuliaStats/HypothesisTests.jl/issues/370). Only the
p-value is affected. `confint` still returns an interval, since
[§6.3](@ref "6.3 Exact index") inverts the lattice distribution whether or not the sample is
tied. Its own cost is bounded separately, as
[§9](@ref "9. The rank tests in this package") sets out; the limits in the table above are
on the p-value alone, and the untied lattice carries none.

The interval follows the type rather than this rule: `Exact*` inverts the lattice
distribution ([§6.3](@ref "6.3 Exact index")), tied sample or not, and `Approximate*` the
normal one ([§6.4](@ref "6.4 Normal-approximation index")), whatever distribution the
p-value beside it was read off.

#### 2.2.1 The lattice distribution

So called because with no ties the statistic lives on the integer lattice: the midranks
are ``1, \dots, n`` whatever the data, and ``W^+`` is a sum over a subset of them.

Under the null the signs of the ``D_i`` are independent, each ``\pm`` with probability
``1/2``, and independent of ``|D|``. With no ties the midranks are the integers
``1, \dots, n``, so

```math
W^+ \;\overset{d}{=}\; \sum_{i=1}^{n} i \, B_i , \qquad B_i \overset{\text{iid}}{\sim} \mathrm{Bernoulli}(1/2) .
```

Let ``c_n(w) = \#\{ S \subseteq \{1,\dots,n\} : \sum_{i \in S} i = w \}``. Then

```math
c_n(w) = c_{n-1}(w) + c_{n-1}(w - n), \qquad
c_0(0) = 1, \quad c_0(w) = 0 \ (w \ne 0), \quad c_n(w) = 0 \ (w < 0),
```

and ``\mathbb{P}(W^+ = w) = c_n(w) / 2^n``. Evaluating the recursion over the whole support is
``O(n^3)`` in time and ``O(n^2)`` in space. Write ``F_n(w) = \mathbb{P}(W^+ \le w)``.

!!! warning "Numerical care"
    ``c_n(w)`` grows like ``2^n``, so accumulating the counts and dividing only at the end
    overflows a fixed-width integer at ``n = 72`` in 64 bits, and does so silently: the
    largest count reaches ``1.05 \times 10^{19}`` against a ceiling of
    ``9.22 \times 10^{18}``. Carrying the normalisation inside the recursion avoids it,
    propagating probabilities in place of counts. Exact arithmetic is the alternative.

    This package takes the first route, through `StatsFuns.signrankcdf`, and
    [§9](@ref "9. The rank tests in this package") records the version from which that holds.

**Implementations.** Here ``F_n`` is `StatsFuns.signrankcdf(n, w)`, with `signrankccdf`
for the upper tail. R computes the same distribution in the C routines `psignrank`,
`dsignrank` and `qsignrank` of its `nmath` library, which accumulate the counts ``c_n(w)``
in `double`: exact while they stay under ``2^{53}`` and approximate beyond, which is a
different way out of the overflow above. The Julia values are tested against R's, as are
the p-values and intervals of [§8](@ref "8. Worked values for the rank tests").

#### 2.2.2 The permutation distribution

So called after the permutation-test literature, whose distributions this is one of:
[Streitberg and Röhmel (1986)](@cite streitberg1986) and `exactRankTests` both carry the
name, the latter in its distribution functions `dperm`, `pperm` and `qperm`. For one sample
the "permutations" are sign flips rather than rearrangements, and the name travels across
anyway. Its defining property is that it is built *conditionally on the midranks the sample
actually produced*, which under ties differ from sample to sample.

With ties the midranks are not ``1, \dots, n``, so the lattice recursion of
[§2.2.1](@ref "2.2.1 The lattice distribution"), which counts subsets of ``\{1, \dots, n\}``, no
longer describes the statistic. What is done instead is to build the null distribution
directly from the sample in hand, in three steps:

1. Rank the observed ``|d_i|``, ties included, giving midranks ``R_1, \dots, R_n``. These
   are now fixed numbers, not ``1, \dots, n``.
2. Under the null each sign is ``\pm`` with probability ``1/2`` independently of ``|D|``,
   so all ``2^n`` ways of attaching signs to those fixed midranks are equally likely. Walk
   through all of them, and for each compute the value ``W^+`` would have taken:
   ``\sum_i \varepsilon_i R_i`` with ``\varepsilon_i \in \{0,1\}`` saying whether
   observation ``i`` was counted as positive.
3. The resulting multiset of ``2^n`` values, each with weight ``2^{-n}``, *is* the null
   distribution. Tail probabilities are counts of it.

That is, the distribution is not looked up or approximated; it is tabulated by brute force
from the observed midranks, which is why ties cost ``O(2^n)`` where their absence costs a
polynomial recursion. Formally,

```math
\mathbb{P}(W^+ \le w) = 2^{-n} \, \#\Bigl\{ \varepsilon \in \{0,1\}^n : \textstyle\sum_i \varepsilon_i R_i \le w \Bigr\} .
```

The proportions of assignments falling at or below and at or above the
observed ``W^+`` are written ``q_{\le}`` and ``q_{\ge}``, and are what
[§2.3](@ref "2.3 p-values") uses.

**Worked, in full.** Take ``d = (1.5,\, -1.5,\, 2.0)``. The absolute values are
``(1.5, 1.5, 2.0)``: the two ``1.5``s share positions ``1`` and ``2`` and so take the
midrank ``1.5`` each, and ``2.0`` takes ``3``. Observed
``W^+ = 1.5 + 3 = 4.5``, the first and third observations being the positive ones. All
``2^3 = 8`` sign assignments to the fixed midranks ``(1.5, 1.5, 3)``:

| signs of ``(d_1, d_2, d_3)`` | midranks counted | ``W^+`` |
|:---|:---|:---|
| ``-\,-\,-`` | none | ``0`` |
| ``+\,-\,-`` | ``1.5`` | ``1.5`` |
| ``-\,+\,-`` | ``1.5`` | ``1.5`` |
| ``-\,-\,+`` | ``3`` | ``3`` |
| ``+\,+\,-`` | ``1.5 + 1.5`` | ``3`` |
| ``+\,-\,+`` | ``1.5 + 3`` | ``4.5`` |
| ``-\,+\,+`` | ``1.5 + 3`` | ``4.5`` |
| ``+\,+\,+`` | ``1.5 + 1.5 + 3`` | ``6`` |

so the null distribution is ``0, 1.5, 3, 4.5, 6`` with probabilities
``\tfrac{1}{8}, \tfrac{2}{8}, \tfrac{2}{8}, \tfrac{2}{8}, \tfrac{1}{8}``. The tie has
collapsed what would have been six distinct rank sums onto five values, which is exactly
how it enters. At the observed ``W^+ = 4.5``, seven of the eight assignments are at or
below it and three are at or above, so ``q_{\le} = 7/8 = 0.875`` and
``q_{\ge} = 3/8 = 0.375``, giving the two-sided ``p = 2 \times 0.375 = 0.75`` by
[§2.3](@ref "2.3 p-values"). Note the distribution is symmetric about
``\tfrac{1}{2}\sum_i R_i = 3``, as it is for any tie pattern: flipping every sign sends
``W^+`` to ``\sum_i R_i - W^+``.

**What it conditions on.** The permutation distribution is a conditional distribution,
and the conditioning is worth stating precisely: it is not the distribution
of ``W^+`` over repeated samples from ``F``. It is the distribution over the ``2^n`` sign
patterns with the observed absolute values, and therefore their midranks, held fixed at
what was seen. Every probability in
[§2.3](@ref "2.3 p-values") is computed in that fixed-``|d|`` distribution.

Two things make this the right object rather than a retreat from one. Under the null the
signs are independent of ``|D|``, so fixing ``|d|`` discards nothing that bears on the
null. And a test whose level is exactly ``\alpha`` for every possible value of ``|d|`` has
level exactly ``\alpha`` when averaged over ``|d|``, which is the unconditional statement:
conditional exactness is the stronger property, not a weaker substitute for it. That is
also why the untied case of [§2.2.1](@ref "2.2.1 The lattice distribution") needs no
such discussion. There the midranks are ``1, \dots, n`` whatever the data, so conditioning
on them fixes nothing, and the two distributions coincide.

Discarding zeros ([§2.1](@ref "2.1 Model, estimand, statistic")) conditions in exactly the
same way, on which observations were non-zero, and is exact for the same reason.

**Implementations.** This package computes this distribution itself, and reaches it on the
routing of [§2.2](@ref "2.2 Null distribution of the signed rank statistic"): every tied
sample it accepts as exact comes here. The tied p-value of
[§8.2](@ref "8.2 One sample, five zeros and ties among the rest") is computed this way.

Nothing in StatsFuns covers the tied case, and base R declines it: `wilcox.test` warns that
it cannot compute an exact p-value with ties, or with zeros, and falls back on its normal
approximation. R's
[`exactRankTests`](https://cran.r-project.org/package=exactRankTests)`::wilcox.exact` does
compute it, and this package's tied p-values are tested against it.

It does so by a better algorithm, worth naming because the ``O(2^n)`` above is not intrinsic
to the problem. ``W^+ = \sum_i \varepsilon_i R_i`` with ``\varepsilon \in \{0,1\}^n`` is a
subset-sum count, so the counts that make up the distribution are the coefficients of a
polynomial: expanding ``\prod_i (1 + z^{2R_i})`` chooses each sign once per factor, and the
coefficient of ``z^{2w}`` is the number of sign patterns with ``W^+ = w``, the whole
distribution after division by ``2^n``. The doubled exponents keep the half-integer
midranks integral. Building
that product by convolving one factor at a time, each step shifting the running counts by
``2R_i`` and adding, gives every coefficient in ``O(n \sum_i R_i) = O(n^3)`` operations
rather than ``O(2^n)``: the **shift algorithm** of
[Streitberg and Röhmel (1986)](@cite streitberg1986). The gain is not marginal. On a tied
sample of ``n = 24`` it returns the same two tails some forty times faster, and at
``n = 60``, where enumerating ``2^{60} \approx 1.2 \times 10^{18}`` sign patterns is out of
the question, it still finishes in milliseconds. This package enumerates, and bounds the
enumeration instead ([§9](@ref "9. The rank tests in this package")); adopting the shift
algorithm would remove that bound rather than raise it, which is proposed as
[issue #370](https://github.com/JuliaStats/HypothesisTests.jl/issues/370).
`wilcox.exact`'s one-sided values agree with [§2.3](@ref "2.3 p-values") exactly; so do its
two-sided ones
here, but for a reason worth recording: `wilcox.exact` doubles nothing, instead summing the
attainable outcomes at least as extreme on the far side, and for the one-sample statistic
the two rules coincide because the one-sample permutation distribution is symmetric
whatever the ties,
the flip of every sign taking ``W^+`` to ``\sum_i R_i - W^+``. The two-sample permutation
distribution has no such symmetry, and there the two rules part:
[§3.2.2](@ref "3.2.2 The permutation distribution").

#### 2.2.3 The normal approximation

``W^+`` is asymptotically normal with the mean and variance of
[§2.1](@ref "2.1 Model, estimand, statistic"). Write

```math
\mu = W^+ - \frac{n(n+1)}{4}, \qquad
\sigma = \sqrt{\frac{n(n+1)(2n+1)}{24} - \frac{T(|d|)}{48}}
```

for the centred statistic and the tie-corrected standard deviation. The variance
correction is exact under the permutation distribution of
[§2.2.2](@ref "2.2.2 The permutation distribution"), not an approximation to it.

### 2.3 p-values

Exact, no ties:

```math
p_{\text{left}} = F_n(W^+), \qquad
p_{\text{right}} = 1 - F_n(W^+ - 1), \qquad
p_{\text{both}} = \min\bigl(1,\, 2\min(p_{\text{left}},\, p_{\text{right}})\bigr) .
```

``F_n`` is symmetric about ``n(n+1)/4``, so the smaller tail is the left one exactly when
``W^+ \le n(n+1)/4``; the two-sided value may equivalently be computed by branching on that
comparison and doubling the selected tail.

The clip is not redundant. Where ``n(n+1)/2`` is even the null mean is itself attainable,
and at ``W^+ = n(n+1)/4`` the two tails are equal, each exceeding ``1/2`` by half the atom
sitting on the mean:

```math
2 F_n\bigl(n(n+1)/4\bigr) = 1 + \mathbb{P}\bigl(W^+ = n(n+1)/4\bigr) > 1 .
```

At ``n = 3`` with ``W^+ = 3`` that doubled tail is ``1.25``. Doubling a discrete tail can
overshoot, and this is where it does.

Exact, ties present: with ``W'`` the statistic recomputed under a sign assignment, and
``q_{\le}`` and ``q_{\ge}`` the proportions of assignments giving ``W' \le W^+`` and
``W' \ge W^+`` respectively,

```math
p_{\text{left}} = q_{\le}, \qquad p_{\text{right}} = q_{\ge}, \qquad
p_{\text{both}} = \min\bigl(1,\, 2\min(q_{\le},\, q_{\ge})\bigr) .
```

Normal approximation, with a continuity correction of ``1/2``:

```math
p_{\text{left}} = \Phi\!\left(\frac{\mu + 1/2}{\sigma}\right), \qquad
p_{\text{right}} = 1 - \Phi\!\left(\frac{\mu - 1/2}{\sigma}\right), \qquad
p_{\text{both}} = 2\left[1 - \Phi\!\left(\frac{\bigl|\mu - \tfrac{1}{2}\operatorname{sign}\mu\bigr|}{\sigma}\right)\right] .
```

Degenerate cases: if ``n = 0``, meaning every difference was zero, all three exact p-values are
``1``. If ``\mu = \sigma = 0``, all three approximate p-values are ``1``.

## 3. The two-sample procedure (Wilcoxon rank sum, Mann-Whitney U)

The two names are one procedure, arrived at independently by
[Wilcoxon (1945)](@cite wilcoxon1945) and [Mann and Whitney (1947)](@cite mann1947).
Wilcoxon tabulated the rank sum ``W = \sum_{i=1}^{n_x} R_i`` itself; ``U`` subtracts its
minimum. Absent ties that minimum is ``n_x(n_x+1)/2``, since the ``n_x`` ranks going to
``x`` are then distinct members of ``\{1, \dots, N\}`` and so total at least
``1 + 2 + \dots + n_x``, attained exactly when every ``x_i`` falls below every ``y_j``.
Under ties it is the same number, which is easier to see from the counting form of
[§3.1](@ref "3.1 Model, estimand, statistic"), where ``U`` counts pairs and so cannot be
negative. Hence ``W = U + n_x(n_x+1)/2``, a shift by a constant fixed by ``n_x`` alone and
never by the data. Being a strictly increasing bijection, it leaves the ordering of outcomes
untouched, so the two statistics give the same tail events, the same null distribution up to
the shift, and the same p-value at every level. Subtracting the minimum is what turns the
rank sum into the pair count of [§1.2](@ref "1.2 A worked ranking"), supported from ``0``
and, divided by ``n_x n_y``, an estimate of
``\mathbb{P}(X > Y) + \tfrac{1}{2} \mathbb{P}(X = Y)``.

The three subsections are those of
[§2](@ref "2. The one-sample procedure (Wilcoxon signed rank)"): the statistic
([§3.1](@ref "3.1 Model, estimand, statistic")), its null distribution
([§3.2](@ref "3.2 Null distribution of the Mann-Whitney statistic")), and the p-values
([§3.3](@ref "3.3 p-values")).

### 3.1 Model, estimand, statistic

**Model.** The two samples are independent random variables ``X_1, \dots, X_{n_x}`` and
``Y_1, \dots, Y_{n_y}``, independent of each other and i.i.d. within each sample. The null
hypothesis is ``F_x = F_y``, the two distributions equal and otherwise unrestricted. What that
requires, and what it does not, is [§3.1.1](@ref "3.1.1 Assumptions").

**Estimand.** Under the **shift model** ``F_x(t) = F_y(t - \Delta)``, the estimand is
``\Delta``. Without a shift model the null tested is still equality of the two
distributions, and what [§5](@ref "5. Point estimation") estimates is the median of
``X - Y`` for independent ``X \sim F_x`` and ``Y \sim F_y``. That quantity is not in general
``\operatorname{median}(F_x) - \operatorname{median}(F_y)``. The direction the test is
sensitive to is a departure from ``\mathbb{P}(X > Y) = 1/2``, which is not the same as that equality
being a null it holds its level against: [§3.1.1](@ref "3.1.1 Assumptions") measures the
difference. Zero observations carry no
special meaning here and are not discarded.

**Hypotheses.** As in [§2.1](@ref "2.1 Model, estimand, statistic"), the null the level is
exact for is a statement about distributions,

```math
H_0 : F_x = F_y ,
```

which under the shift model ``F_x(t) = F_y(t - \Delta)`` reads ``\Delta = 0``, with

```math
H_1 : \Delta \ne 0 \quad (\texttt{tail = :both}), \qquad
H_1 : \Delta < 0 \quad (\texttt{tail = :left}), \qquad
H_1 : \Delta > 0 \quad (\texttt{tail = :right}) ,
```

small ``U`` being the evidence that ``x`` sits below ``y``. The one-sided intervals of
[§6.5](@ref "6.5 One-sided intervals") take the same names the same way. Testing
``\Delta = \Delta_0`` means running the procedure on ``x_i - \Delta_0`` against ``y``, which
is the ``U(\Delta)`` of [§4.2](@ref "4.2 The counting identity").

**Statistic.** Rank the pooled sample of size ``N``. With ``R_i`` the midrank of ``x_i``,

```math
U = \sum_{i=1}^{n_x} R_i - \frac{n_x(n_x+1)}{2}
  = \#\{(i,j) : x_i > y_j\} + \tfrac{1}{2} \#\{(i,j) : x_i = y_j\} .
```

The two forms are algebraically identical, and the subtracted ``n_x(n_x+1)/2`` is what
makes them so. Split the pooled midrank of ``x_i`` according to which sample each
counted value came from:

```math
R_i = \left[ \#\{j : x_j < x_i\} + \frac{1 + \#\{j : x_j = x_i\}}{2} \right]
    + \#\{j : y_j < x_i\} + \tfrac{1}{2} \#\{j : y_j = x_i\} .
```

The bracket is the midrank of ``x_i`` within ``x`` alone, so summing it over ``i`` gives
``n_x(n_x+1)/2`` however ``x`` is tied
([§1.1](@ref "1.1 Mathematical observations")), which is precisely the term subtracted.
What survives is the pair count. The first form is how the statistic is computed, from one
sort of the pooled sample; the second is the definition the inversion of
[§6.2](@ref "6.2 Inversion") works with. The support is
``\{0, \tfrac{1}{2}, \dots, n_x n_y\}``, integer-valued absent ties, symmetric under the
null about ``n_x n_y / 2``.

**Which sample is called ``x``.** The labelling fixes the sign of what is estimated and
which tail is which; nothing that is tested depends on it. Write ``U_x`` for the statistic
above and ``U_y`` for the one obtained by
exchanging the roles of the two samples. Every pair ``(i,j)`` is counted once between them,
as ``x_i > y_j``, as ``y_j > x_i``, or as a tie splitting half to each, so

```math
U_x + U_y = n_x n_y
```

identically, ties included and whatever the data. Exchanging the samples therefore reflects
the statistic through ``n_x n_y / 2``, which is exactly the centre of symmetry of its null
distribution. So the two-sided p-value is unchanged, and the one-sided ones swap: the left
tail of ``U_x`` is the right tail of ``U_y``. On the two-sample example of
[§1.2](@ref "1.2 A worked ranking"), ``U_x = 6.5`` and ``U_y = 13.5``, both giving
``p = 0.44444`` two-sided, and ``p = 0.22222`` for the left tail of the first and the right
tail of the second. Interval and estimate follow the sign of the pairwise estimates of
[§4](@ref "4. Pairwise estimates"): exchanging the samples negates every ``D_{ij}``, and so
negates and reverses both, taking the interval of
[§8.4](@ref "8.4 Two samples, ties") from ``(-14, -1)`` to ``(1, 14)``.

That the two are interchangeable is what lets an implementation work with whichever is
convenient, and this package uses the freedom: the tied enumeration of
[§3.2.2](@ref "3.2.2 The permutation distribution") always enumerates the *smaller*
sample, since there are fewer subsets to visit, and swaps the two tails afterwards if that
was ``y``.

Under the null the pooled midranks are fixed and the ``n_x`` of them falling to ``x`` are
a simple random sample without replacement, so ``\sum_{i} R_i`` over that sample has
variance ``\frac{n_x n_y}{N(N-1)}\sum_i (R_i - \bar R)^2``. Substituting
[§1.1](@ref "1.1 Mathematical observations"),

```math
\mathbb{E}[U] = \frac{n_x n_y}{2}, \qquad
\operatorname{Var}(U) = \frac{n_x n_y}{12}\left(N + 1 - \frac{T([x; y])}{N(N-1)}\right) .
```

#### 3.1.1 Assumptions

Equality is what makes the test exact, because it makes the ``N`` observations
exchangeable: every assignment of the pooled midranks to the two samples is then equally
likely, which is the whole of [§3.2.1](@ref "3.2.1 The lattice distribution") and
[§3.2.2](@ref "3.2.2 The permutation distribution"). As in
[§2.1.1](@ref "2.1.1 Assumptions"), continuity of ``F_x`` and ``F_y`` is a
convenience: it makes ties events of probability zero, so the lattice distribution of
[§3.2.1](@ref "3.2.1 The lattice distribution") applies. Discrete or mixed
distributions keep exactness through the permutation enumeration of
[§3.2.2](@ref "3.2.2 The permutation distribution").

Equality cannot be weakened to ``\mathbb{P}(X > Y) = 1/2``. That weaker statement leaves ``F_x`` and
``F_y`` free to differ in spread, and then the pooled observations are no longer
exchangeable, the null variance above is no longer the variance of ``U``, and the
test does not hold its level. It fails in both directions, according to which sample is
the larger: for ``X \sim \mathcal{N}(0, 1)`` against ``Y \sim \mathcal{N}(0, 9)``, where
``\mathbb{P}(X > Y) = 1/2`` holds exactly by symmetry, the nominal ``0.05`` two-sided test has size
about ``0.13`` at ``(n_x, n_y) = (30, 10)`` and about ``0.016`` at ``(10, 30)``. Under
``F_x = F_y`` both come to ``0.05``. This is the two-sample counterpart of the median
against pseudomedian trap in [§2.1.1](@ref "2.1.1 Assumptions"), and it is why
that section insists on symmetry and this one on equality.

### 3.2 Null distribution of the Mann-Whitney statistic

One exact null distribution and the same three forms as in
[§2.2](@ref "2.2 Null distribution of the signed rank statistic"), carrying the same names,
with the enumeration now over rank assignments rather than sign patterns. Write
``B = \binom{N}{\min(n_x, n_y)}`` for the number of those assignments:

| name | what it is | available for | cost |
|:---|:---|:---|:---|
| **lattice** ([§3.2.1](@ref "3.2.1 The lattice distribution")) | the exact distribution in its untied form, tabulated by a recursion | untied pooled samples only, ``T([x; y]) = 0`` | ``O\bigl((n_x n_y)^2\bigr)`` |
| **permutation** ([§3.2.2](@ref "3.2.2 The permutation distribution")) | the exact distribution in general, enumerated from the sample | any sample, and the only exact form once ``T([x; y]) > 0`` | ``O(B)`` |
| **normal** ([§3.2.3](@ref "3.2.3 The normal approximation")) | an approximation to it, from the moments above | any sample | ``O(N)`` |

**Choosing between them.** Exactly as in
[§2.2](@ref "2.2 Null distribution of the signed rank statistic"): the tie pattern decides
whether the lattice form is available, and exact against normal is a choice about cost, so
an untied sample too large for the recursion still goes through the normal. Here too
lattice and permutation are two computations of the one exact distribution, coinciding
whenever the pooled midranks come out as ``1, \dots, N``.

**In this package.** As in
[§2.2](@ref "2.2 Null distribution of the signed rank statistic"), two types implement the
procedure, `ExactMannWhitneyUTest` and `ApproximateMannWhitneyUTest`, and `MannWhitneyUTest`
is a function returning one of them rather than a type of its own. The tie pattern and the
pooled size are the only properties of the data that enter:

| call | pooled sample untied | pooled sample tied |
|:---|:---|:---|
| `MannWhitneyUTest(x, y)` | lattice for ``N \le 50``, normal above | permutation for ``N \le 10``, normal above |
| `ExactMannWhitneyUTest(x, y)` | lattice, no bound | permutation, up to ``B = 2^{30}`` |
| `ApproximateMannWhitneyUTest(x, y)` | normal | normal |

Again the automatic rule turns to the normal approximation once the sample outgrows its exact
route, and ties lower the size at which that happens, here from ``N \le 50`` to ``N \le 10``.
``N`` counts every observation, since nothing is discarded
([§3.1](@ref "3.1 Model, estimand, statistic")).

The two tied thresholds, ``n \le 15`` there and ``N \le 10`` here, are long-standing choices
of this package rather than two readings of one cost target, and they are not equally strict.
At the signed rank limit the enumeration visits ``2^{15} = 32\,768`` sign patterns; at the
Mann-Whitney limit it visits at most ``\binom{10}{5} = 252`` assignments. The two-sample rule
therefore turns to the approximation a good deal earlier than its own cost requires.

`method` overrides the rule as it does in
[§2.2](@ref "2.2 Null distribution of the signed rank statistic"), the callable here being
passed `(; nx, ny, ties, tie_adjustment)`. Forcing the exact route is declined the same way,
here when ``B`` exceeds ``2^{30} = 1\,073\,741\,824`` assignments: `pvalue` throws
[`ComputationTooLarge`](@ref HypothesisTests.ComputationTooLarge) naming ``n_x`` and ``n_y``.
This bound does not reach `confint`, which inverts the untied lattice whichever route the
p-value took, but the exact interval carries a bound of its own: see
[§9](@ref "9. The rank tests in this package").

#### 3.2.1 The lattice distribution

Under the null all ``\binom{N}{n_x}`` assignments of pooled ranks to the two samples are
equally likely. Let ``c_{n_x, n_y}(u)`` count those giving ``U = u``. Conditioning on
whether the largest pooled observation belongs to ``x`` or to ``y``,

```math
c_{n_x, n_y}(u) = c_{n_x - 1,\, n_y}(u - n_y) + c_{n_x,\, n_y - 1}(u) ,
```

with ``c_{n_x, 0}(0) = c_{0, n_y}(0) = 1``, and ``c(u) = 0`` for ``u < 0`` or
``u > n_x n_y``. Then ``\mathbb{P}(U = u) = c_{n_x,n_y}(u) / \binom{N}{n_x}``. Evaluating over the
support is ``O\bigl((n_x n_y)^2\bigr)`` in time and ``O(n_x n_y)`` in space. Write
``G_{n_x,n_y}(u) = \mathbb{P}(U \le u)``.

The numerical caveat of [§2.2.1](@ref "2.2.1 The lattice distribution") applies with
more force: the normalising constant ``\binom{N}{n_x}`` exceeds ``2^{63}`` for balanced
samples from ``n_x = n_y = 34``, where ``\binom{68}{34} \approx 2.85 \times 10^{19}``.
That bounds the counts themselves; an implementation whose intermediate terms exceed the
counts overflows earlier still, which is why carrying the normalisation inside the recursion
matters more here than it does for one sample.

**Implementations.** Here ``G_{n_x,n_y}`` is `StatsFuns.wilcoxcdf(nx, ny, u)`, with
`wilcoxccdf` for the upper tail. It computes this distribution, but not by the recursion
above: it uses the recurrence of [Löffler (1983)](@cite loeffler1983), a convolution over
divisor sums that reaches the same probabilities with fewer allocations. The recursion
stated here is the one R uses, in the C routines `pwilcox`, `dwilcox` and `qwilcox` of its
`nmath` library, which memoise the counts in `double`. Two algorithms, one distribution;
that the two agree is worth checking, and they do, to floating-point precision at every
attainable ``u`` for ``n_x, n_y \le 7``.

Ties are handled as in [§2.2.2](@ref "2.2.2 The permutation distribution"), by this package's own
enumeration, on the routing of
[§3.2](@ref "3.2 Null distribution of the Mann-Whitney statistic"). Neither StatsFuns nor
base R offers it; `exactRankTests` does, as
[§2.2.2](@ref "2.2.2 The permutation distribution") records.

#### 3.2.2 The permutation distribution

The recipe of [§2.2.2](@ref "2.2.2 The permutation distribution") carries over with the sign
patterns replaced by group assignments. Rank the pooled sample, ties included, fixing the
``N`` midranks. Under the null the two samples are exchangeable, so every way of dealing
``\min(n_x, n_y)`` of those fixed midranks to the smaller sample is equally likely; walk
through all ``\binom{N}{\min(n_x, n_y)}`` of them, compute the ``U`` each would give, and
the resulting multiset, each member weighted equally, is the null distribution. The same
three p-value formulas then follow.

On ``x = (1, 2)`` against ``y = (2, 3, 4)`` the pooled midranks are
``(1,\, 2.5,\, 2.5,\, 4,\, 5)``, the two ``2``s sharing positions ``2`` and ``3``.
Dealing two of them to ``x`` in all ``\binom{5}{2} = 10`` ways gives

```math
U \in \{0.5,\, 2,\, 3,\, 3.5,\, 4.5,\, 6\}
\quad\text{with counts}\quad 2,\, 2,\, 1,\, 2,\, 2,\, 1 \ \text{out of } 10 ,
```

and the observed ``U = 0.5``.

That distribution is worth a second look, because unlike its one-sample counterpart it is
in general *not* symmetric. Here the centre would be ``n_x n_y / 2 = 3``, yet ``0.5`` carries
weight ``2/10`` while its reflection ``5.5`` carries none. Symmetry would need the tie
pattern to read the same from both ends, which nothing guarantees.

The one-sided p-values are unaffected by that. The two-sided one is where
conventions part, since doubling the smaller tail and summing the far tail now disagree:
on ``x = 1, \dots, 10`` against ``y = 2, 4, \dots, 24``, the doubling of
[§3.3](@ref "3.3 p-values") gives ``0.0120035`` where `exactRankTests::wilcox.exact`,
which sums, gives ``0.0118890``, the one-sided values agreeing exactly. Neither is an error;
they answer differently the question of what counts as at least as extreme on the other side
of an asymmetric distribution. This specification doubles, which is also what the untied
routes of both procedures do.

#### 3.2.3 The normal approximation

``U`` is asymptotically normal with the mean and variance of
[§3.1](@ref "3.1 Model, estimand, statistic"). Write

```math
\mu = U - \frac{n_x n_y}{2}, \qquad
\sigma = \sqrt{\frac{n_x n_y}{12}\left(N + 1 - \frac{T([x; y])}{N(N-1)}\right)}
```

for the centred statistic and the tie-corrected standard deviation. As in
[§2.2.3](@ref "2.2.3 The normal approximation"), the tie correction is exact under the permutation
distribution of [§3.2.2](@ref "3.2.2 The permutation distribution") rather than an
approximation to it; what is approximate is only the normal shape.

### 3.3 p-values

Exact, no ties. Write ``G = G_{n_x,n_y}`` for the null cdf of
[§3.2.1](@ref "3.2.1 The lattice distribution") and ``U'`` for a draw from it, so that
``G(u) = \mathbb{P}(U' \le u)``; the sample sizes are fixed throughout, so carrying them as
subscripts adds nothing here. Then

```math
p_{\text{left}} = G(U), \qquad
p_{\text{right}} = 1 - G(U - 1), \qquad
p_{\text{both}} = \min\bigl(1,\, 2\, G(\min(U,\, n_x n_y - U))\bigr) .
```

The fold to the lower tail is exact, not an approximation: by the null symmetry of
[§3.1](@ref "3.1 Model, estimand, statistic"),
``G(n_x n_y - U) = \mathbb{P}(U' \ge U)``, so ``G(\min(U, n_x n_y - U))`` is the smaller of the two
tails wherever ``U`` sits, the centre included. The clip is needed for the reason it is
needed in [§2.3](@ref "2.3 p-values"): what can exceed ``1`` is the doubling, when ``U``
lands on ``n_x n_y / 2`` and that value carries an atom. At ``n_x = n_y = 2`` with
``U = 2`` the doubled tail is ``4/3``.

Exact under ties, and the normal approximation: exactly as in [§2.3](@ref "2.3 p-values"),
with ``q_{\le}, q_{\ge}`` from
[§3.2.2](@ref "3.2.2 The permutation distribution") and
``\mu, \sigma`` from [§3.2.3](@ref "3.2.3 The normal approximation").

Degenerate cases: if ``n_x = 0`` or ``n_y = 0`` there is no pair to compare, ``U = 0`` is
the only attainable value, the null distribution is a point mass there, and all three exact
p-values are ``1``. If ``\mu = \sigma = 0``, all three approximate p-values are ``1``.

## 4. Pairwise estimates

The two procedures differ in one place only, and this section is where that difference is
isolated. Each forms a set of ``m`` numbers on the scale of the data, one for each pair of
observations its statistic compares. Each of those numbers is, on its own, a crude estimate
of the quantity the procedure is after: ``\theta`` from one pair of differences,
``\Delta`` from one ``x`` against one ``y``. This specification therefore calls them the
**pairwise estimates**. What differs between the procedures is which pairs there are, and
what number each pair contributes:

  - the one-sample statistic compares every retained observation with every other and with
    itself, giving the ``m = n(n+1)/2`` pairs ``i \le j`` of a single sample, each
    contributing the average ``(d_i + d_j)/2``;
  - the two-sample statistic compares every ``x`` with every ``y``, and never an ``x`` with
    an ``x``, giving the ``m = n_x n_y`` cross pairs, each contributing the difference
    ``x_i - y_j``.

[§4.1](@ref "4.1 Definitions") states those two definitions, and they are the last thing
stated twice for a substantive reason. Everything after them is a statement about a set of
``m`` numbers: sort them as ``V_{(1)} \le \dots \le V_{(m)}``, and the estimate of
[§5](@ref "5. Point estimation") is their sample median while the interval endpoints of
[§6](@ref "6. Interval estimation") are two of their order statistics. Pooling ``m`` crude
estimates by taking their median is the whole of what the procedure does with them.

Note which median that is. The median of the ``m`` pairwise estimates is a number computed
from the data, and it is the estimator; the pseudomedian of
[§2.1](@ref "2.1 Model, estimand, statistic") and the shift of
[§3.1](@ref "3.1 Model, estimand, statistic") are the population quantities being
estimated. They are not the same object, and neither is the sample median of the data
itself. [§5](@ref "5. Point estimation") keeps the three apart.

The common treatment rests on the identity of
[§4.2](@ref "4.2 The counting identity"): for both procedures the statistic, recomputed
against a hypothesised location, counts how many pairwise estimates lie above it. Inverting
either test is therefore counting them, and that is the only property
[§5](@ref "5. Point estimation") and [§6](@ref "6. Interval estimation") use.

The notation below nevertheless stays doubled, because the two statistics keep their own
names: results are stated for the one-sample case, in ``W^+`` and ``\theta``, and
[§6.2](@ref "6.2 Inversion") gives the substitution (``U``, ``D``, ``\Delta``) that carries
each of them to the two-sample case. No step of any argument changes under that
substitution. Reading the one-sample line and applying the substitution is enough; the two
are not being developed independently.

The term **pairwise estimates** is a convenience of this document rather than established
usage. The one-sample ones are normally called the Walsh averages and the two-sample ones
simply the differences, and no standard term covers both. Where the literature does name the
common idea it speaks of *elementary estimates*, the same thing: one estimate per pair,
combined by taking a median.

### 4.1 Definitions

**One sample: the Walsh averages.** The ``m = n(n+1)/2`` values

```math
A_{ij} = \frac{d_i + d_j}{2}, \qquad 1 \le i \le j \le n ,
```

over the ``n`` retained observations. The diagonal ``i = j`` is included, so each ``d_i``
is itself a member.

**Two samples: the cross-group differences.** The ``m = n_x n_y`` values

```math
D_{ij} = x_i - y_j , \qquad 1 \le i \le n_x, \ 1 \le j \le n_y .
```

``D_{ij}`` is the traditional letter for these, and its two indices are what set it apart
from the random variable ``D_i`` of §1: a cross
difference is a number computed from the two samples, not an observation regarded as
random.

**Cost.** Both sets are quadratic in the sample size, and this specification forms them
explicitly: ``m`` numbers, sorted. That is the binding constraint at scale, and it is
avoidable, because the estimate of [§5](@ref "5. Point estimation") and the interval
endpoints of [§6](@ref "6. Interval estimation") are order statistics. A selection
algorithm returns them without materialising the set, exploiting the fact that the Walsh
averages and the cross differences each form a sorted matrix. An implementation that
materialises instead should bound the size it will accept.

### 4.2 The counting identity

Let ``W^+(\theta)`` denote the statistic of [§2.1](@ref "2.1 Model, estimand, statistic")
recomputed on ``d_1 - \theta, \dots, d_n - \theta``, and ``U(\Delta)`` the statistic of
[§3.1](@ref "3.1 Model, estimand, statistic") recomputed on
``x_1 - \Delta, \dots, x_{n_x} - \Delta`` against ``y``. Then

```math
W^+(\theta) = \#\{(i,j) : i \le j, \ A_{ij} > \theta\} + \tfrac{1}{2}\,\#\{(i,j) : i \le j, \ A_{ij} = \theta\} ,
```
```math
U(\Delta) = \#\{(i,j) : D_{ij} > \Delta\} + \tfrac{1}{2}\,\#\{(i,j) : D_{ij} = \Delta\} .
```

The half-count is the midrank convention, and it vanishes whenever ``\theta`` is not
itself one of the pairwise estimates, which is almost surely the case under continuous
``F`` and is the only case [§6.2](@ref "6.2 Inversion") needs. Both functions are
non-increasing in their argument and change value only at pairwise estimates.

Where ``\theta`` *is* one of them, two cases separate. If ``\theta = (d_i + d_j)/2`` for
``i \ne j``, the identity still describes what the procedure computes on ``d - \theta``: that
sample has ``\lvert d_i - \theta \rvert = \lvert d_j - \theta \rvert`` with opposite signs, and
the midrank the pair shares splits exactly as the half-count does. If ``\theta = d_i``, it
does not: ``d_i - \theta`` is zero, so [§2.1](@ref "2.1 Model, estimand, statistic") discards
it and the remaining ``n-1`` observations are re-ranked. On ``d = (1, 2, 4)`` at
``\theta = 2`` the identity gives ``3.5`` and the procedure gives ``2``.
[§6.2.1](@ref "6.2.1 What happens at the endpoints") is where this matters.

*Proof of the first, for ``\theta`` not itself a pairwise estimate and ``|d|`` untied.*
``R_i = \#\{j : |d_j| \le |d_i|\}``, so
``W^+ = \sum_{i : d_i > 0} \#\{j : |d_j| \le |d_i|\}``. A pair ``\{i, j\}`` is counted
exactly when the one with the larger absolute value is positive, which is exactly when
``d_i + d_j > 0``; the diagonal term ``\{i,i\}`` is counted exactly when ``d_i > 0``.
Applying this to ``d - \theta`` gives the statement. The second is immediate from the
counting form of ``U`` in [§3.1](@ref "3.1 Model, estimand, statistic"). ∎

Everything in [§5](@ref "5. Point estimation") and [§6](@ref "6. Interval estimation") is
a consequence of this identity. It is the reason the estimator and the interval are
functions of the pairwise estimates and not of the ranks.

## 5. Point estimation

The **Hodges-Lehmann estimator** [hodges1963](@cite) is the sample median of the pairwise
estimates of [§4.1](@ref "4.1 Definitions"), which gives one estimator per procedure.

**One sample.** The median of the ``n(n+1)/2`` Walsh averages,

```math
\hat\theta = \operatorname{median}\{ A_{ij} : 1 \le i \le j \le n \} ,
```

a consistent estimator of the pseudomedian of
[§2.1](@ref "2.1 Model, estimand, statistic"). It is not the sample median of the ``d_i``,
which estimates the median of ``F``, a different functional except under symmetry.

**Two samples.** The median of the ``n_x n_y`` cross-group differences,

```math
\hat\Delta = \operatorname{median}\{ D_{ij} : 1 \le i \le n_x,\ 1 \le j \le n_y \} ,
```

a consistent estimator of the shift of [§3.1](@ref "3.1 Model, estimand, statistic"), or
without a shift model of the median of ``X - Y``. It is not
``\operatorname{median}(x) - \operatorname{median}(y)``, which is a different quantity
again.

Both statements about what these are *not* matter in practice, because the gap shows up on
ordinary samples. On the one-sample data of
[§8.1](@ref "8.1 One sample, no ties and no zeros"), ``\hat\theta = 9.675`` against a sample
median of ``10.1``. On a tied nine-against-nine
two-sample set, ``\hat\Delta = 0.56`` against a difference of sample medians of ``0.62``.
Exact symmetry makes each pair agree, but the converse fails: ``d = (1, 3, 3, 8)`` is
symmetric about nothing, and its ten Walsh averages
``1, 2, 2, 3, 3, 3, 4.5, 5.5, 5.5, 8`` have median ``3``, exactly the sample median.
Agreement on a given sample is therefore no evidence of symmetry.

For either, by [§4.2](@ref "4.2 The counting identity"), the estimate is the value at which
the statistic sits closest to its null mean ``m/2``, which is the sense in which it is the
estimator the test induces. For even ``m`` the median is taken as the mean of
``V_{(m/2)}`` and ``V_{(m/2+1)}``, so the estimate need not itself be one of
them.

**Changing scale.** Shifting and stretching the data does the same to the estimate: from
``a d + b`` with ``a > 0`` the estimator returns ``a \hat\theta + b``. Anything more than
that changes the answer. Take logs, estimate, exponentiate, and what comes back is a ratio,
a perfectly good estimate of a different quantity, but not the number the procedure returns
on the untransformed data.

The p-value behaves differently, which is easy to conflate. ``U`` depends on the pooled
data only through its ranks, so any increasing transformation of a two-sample data set
leaves the two-sided p-value bit-for-bit unchanged while moving the estimate. Choosing the
scale is therefore part of specifying the analysis, and it is a choice about the estimate,
not about the test.

## 6. Interval estimation

### 6.1 Form

For an integer ``k \in \{0, 1, \dots, \lceil m/2 \rceil - 1\}``, the two-sided interval is

```math
\bigl(\, V_{(k+1)}, \ V_{(m-k)} \,\bigr) ,
```

a pair of order statistics of the pairwise estimates of [§4](@ref "4. Pairwise estimates").
The range of ``k`` stops where the pair would cross: at ``k = m/2``, reachable only for
even ``m``, the left index would exceed the right. Equivalently ``(V_{(C_\alpha)}, V_{(m+1-C_\alpha)})`` with ``C_\alpha = k+1``, which is
the form used in most of the literature. Only the choice of ``k`` distinguishes the exact
construction ([§6.3](@ref "6.3 Exact index")) from the approximate one
([§6.4](@ref "6.4 Normal-approximation index")).

These intervals carry names. The two-sample one, a pair of order statistics of the
``n_x n_y`` cross-group differences, is the **Moses interval**, after the chapter
L. E. Moses contributed to Walker and Lev's *Statistical Inference* (1953); the one-sample
one, read off the Walsh averages, is usually credited to Tukey. Both are also called
distribution-free, or Hodges-Lehmann, confidence intervals, the latter because
[§5](@ref "5. Point estimation") sits inside them by construction. This page treats them as
one object because [§6.2](@ref "6.2 Inversion") derives both from the same counting
identity. [Hollander and Wolfe (1973)](@cite hollander1973) tabulates both, at pages 27–33
and 68–75.

### 6.2 Inversion

Take the one-sample case; the two-sample case is identical with ``U``, ``D`` and
``\Delta`` throughout. By [§4.2](@ref "4.2 The counting identity"), taking ``\theta`` not
itself a pairwise estimate so that the half-count vanishes, ``W^+(\theta) = \#\{A_{ij} > \theta\}``,
so ``\#\{A_{ij} \le \theta\} = m - W^+(\theta)`` and

```math
W^+(\theta) \ge k+1 \iff \theta < A_{(m-k)} , \qquad
W^+(\theta) \le m-k-1 \iff \theta \ge A_{(k+1)} .
```

The two-sided test with rejection region ``\{W^+ \le k\} \cup \{W^+ \ge m-k\}`` therefore
fails to reject exactly on ``\bigl[A_{(k+1)},\, A_{(m-k)}\bigr)``, for every ``\theta`` the
derivation covers, which is every ``\theta`` that is not itself a pairwise estimate. The
endpoints are the two values it does not cover, and
[§6.2.1](@ref "6.2.1 What happens at the endpoints") takes them up: under continuous ``F``
they are hit with probability zero and the question is idle, but that is a statement about
continuous ``F`` and not about the construction. Since the
null distribution of ``W^+`` is symmetric about ``m/2``, the two rejection tails have
equal probability and, written in the case-free ``V_{(\cdot)}`` of
[§4](@ref "4. Pairwise estimates") since nothing one-sample-specific remains,

```math
\mathbb{P}\bigl(\theta \in (V_{(k+1)}, V_{(m-k)})\bigr) = 1 - 2\,\mathbb{P}(W^+ \le k) ,
```

decreasing in ``k``: larger ``k`` gives a narrower interval and less coverage. The
attainable coverages form a finite set, and no construction of this form can achieve a
value between two of them.

**That is an equality, and only under this section's assumptions.** It is exact, not a
bound, when ``F`` is continuous and symmetric about ``\theta``: continuity is what makes
``\theta`` almost surely not a pairwise estimate and the ``|d_i - \theta|`` almost surely untied, and
symmetry is what gives ``W^+(\theta)`` the null distribution whose tails appear on the
right. Drop either and the equality goes:

  - ties or zeros, which is to say discrete data, split the answer by whether the
    endpoints count, a distinction continuity had made idle
    ([§6.2.1](@ref "6.2.1 What happens at the endpoints")). With the endpoints included
    the interval is conservative: on ``d`` uniform on ``\{-3, \dots, 3\}`` at ``n = 15``,
    where essentially every sample is tied, it covers about ``0.982`` against the ``0.95209`` this
    formula gives. Excluding them, as the open interval this section derives does,
    coverage falls to about ``0.888``, *below* nominal: on data this coarse an endpoint
    lands exactly on ``\theta`` in about ``9\%`` of samples, and an open interval
    excludes its own endpoints. The entire gap between the two figures is that boundary
    event. [§6.2.1](@ref "6.2.1 What happens at the endpoints") and
    [§6.6](@ref "6.6 Zeros, ties, and degeneracy") return to this.
  - asymmetric ``F`` breaks the equality by less, and endpoints play no part. For
    ``\mathrm{Exponential}(1)`` at ``n = 15``, coverage of the pseudomedian is about
    ``0.947``.

The figures are Monte Carlo over ``200\,000`` samples through this package's `confint`;
the same experiment on a continuous symmetric ``F`` returns ``0.95227``, the equality's
``0.95209`` to Monte Carlo accuracy, with open and closed identical.

#### 6.2.1 What happens at the endpoints

The derivation above holds for ``\theta`` that is not one of the pairwise estimates, and the
two endpoints are pairwise estimates by construction. Whether they belong to the interval is
therefore not settled by it, and taking the closure would assert more than the derivation
supports.

At such a ``\theta`` the sample actually tested, ``d - \theta``, is not the untied sample the
interval was read off. As [§4.2](@ref "4.2 The counting identity") sets out,
``\theta = (d_i + d_j)/2`` with ``i \ne j`` leaves two shifted observations equal in absolute
value and opposite in sign, so the p-value there comes from the permutation distribution
of [§2.2.2](@ref "2.2.2 The permutation distribution"); and ``\theta = d_i`` leaves one of them zero,
so it comes from a sample one observation shorter. Neither is the distribution whose quantile
fixed ``k``, so neither need agree with the interval about that point.

They do not on the sample of [§8.1](@ref "8.1 One sample, no ties and no zeros"). Its exact
``0.95`` interval is ``(3.3, 15.5)``, and testing ``\theta`` at each endpoint gives
``p = 0.0493`` at ``3.3`` and ``p = 0.0496`` at ``15.5``, both rejections at ``0.05``, while
just inside either endpoint ``p = 0.0554``. Replicating the ``0.0493`` needs the endpoint at
full precision: the Walsh average printed ``3.3`` is ``3.3000000000000003``, and at a literal
``3.3`` the shifted values that should tie in absolute value miss each other in floating
point, giving ``0.0479`` from an untied sample instead. There the interval behaves as the open one, which
is how [§6.1](@ref "6.1 Form") writes it.

**When this matters.** Under continuous ``F`` it does not: the endpoints are pairwise
estimates of the observed sample, and ``\theta`` coincides with one of them with probability
zero, so the coverage statement above is unaffected. It is discrete and coarsely recorded
data that make the case real, since there a hypothesised value routinely coincides with an
observed one: integer scores and counts, measurements rounded to a grid, and times recorded
on a coarse schedule, such as a time to peak read off scheduled sampling points.

On such data the endpoint question is not a nicety but the dominant term. Measured through
this package's `confint` at ``n = 15`` and nominal ``0.95``, over ``200\,000`` samples per
row:

| ``F`` | endpoints included | endpoints excluded | ``\mathbb{P}(\text{an endpoint} = \theta)`` |
|:---|:---|:---|:---|
| uniform on ``\{-3, \dots, 3\}`` | ``0.98174`` | ``0.88818`` | ``0.094`` |
| uniform on ``\pm\{1, 2, 3\}`` | ``0.98324`` | ``0.87114`` | ``0.112`` |
| standard normal | ``0.95227`` | ``0.95227`` | ``0`` |

Two effects separate. Inverting the untied null distribution for tied data is worth about
``+0.03`` of coverage, which is the conservatism
[§6.6](@ref "6.6 Zeros, ties, and degeneracy") describes and the first column shows.
Excluding the endpoints is worth about ``-0.09`` to ``-0.11``, three times as large and in
the other direction, which is the whole distance between the two coverage columns. The
classical conservative guarantee is a statement about the closed interval,
``\mathbb{P}(V_{(k+1)} \le \theta \le V_{(m-k)}) \ge 1 - \alpha``; a reader of the returned pair
who treats the endpoints as excluded holds an interval that on data like this covers less
than it claims. On coarse data, read the endpoints as included.

### 6.3 Exact index

```math
k = \max\bigl\{\, j \in \{0,\dots,\lceil m/2 \rceil - 1\} \;:\; \mathbb{P}(W \le j) < \alpha/2 \,\bigr\} ,
```

taken as ``0`` when no such ``j`` exists, with ``\mathbb{P}(W \le \cdot)`` the exact null CDF of
[§2.2.1](@ref "2.2.1 The lattice distribution") or
[§3.2.1](@ref "3.2.1 The lattice distribution").

By [§6.2](@ref "6.2 Inversion") the attained coverage is then
``1 - 2\mathbb{P}(W \le k) > 1 - \alpha`` strictly, and the next narrower interval, at ``k+1``,
attains at most ``1-\alpha``. This ``k`` therefore gives the narrowest interval of this
form whose coverage still reaches the nominal level. The excess over ``1-\alpha`` is the
discreteness of the null distribution and is not removable.

Equivalently ``C_\alpha = k+1 = \min\{j : \mathbb{P}(W \le j) \ge \alpha/2\}``, the ``\alpha/2``
quantile of the null distribution, whenever that minimum is at least ``1``. Where it is
``0``, which is the degenerate case of
[§6.6](@ref "6.6 Zeros, ties, and degeneracy"), the two disagree: the convention above
gives ``k = 0`` and so ``C_\alpha = 1``, which is the widest interval this form admits,
whereas the quantile would index outside the pairwise estimates.

!!! note "This is a choice"
    The alternative is to take the attainable coverage *closest* to ``1-\alpha``, which
    may fall below it. That is not a conservative interval, and it is not what is
    specified here. The two rules coincide whenever the attainable coverage immediately
    below nominal is further from it than the one immediately above.

``\mathbb{P}(W \le j)`` is monotone in ``j``, so the condition holds on an initial segment of
``\{0, \dots, \lceil m/2 \rceil - 1\}`` and ``k`` is its last member: a binary search finds it
without evaluating the CDF at every index, which matters because each evaluation runs a
lattice recursion of its own.

### 6.4 Normal-approximation index

The target is the exact critical value ``C_\alpha = \min\{j : \mathbb{P}(W \le j) \ge \alpha/2\}``.
The statistic is supported on a unit lattice, so with ``\mu_0`` the null **mean**, not
the centred statistic of [§2.2.3](@ref "2.2.3 The normal approximation"),

```math
\mathbb{P}(W \le j) \approx \Phi\!\left(\frac{j + 1/2 - \mu_0}{\sigma}\right) .
```

Setting this to ``\alpha/2`` and solving for ``j``,

```math
C_\alpha = \bigl\lceil\, \mu_0 - z_{1-\alpha/2}\,\sigma - \tfrac{1}{2} \,\bigr\rceil ,
\qquad k = C_\alpha - 1 ,
```

clamped to ``\{0,\dots,\lceil m/2 \rceil - 1\}``. In both procedures ``\mu_0 = m/2``: for
one sample ``n(n+1)/4 = m/2``, for two ``n_x n_y / 2 = m/2``.

``\sigma`` is the tie-corrected standard deviation of
[§2.2.3](@ref "2.2.3 The normal approximation") or [§3.2.3](@ref "3.2.3 The normal approximation"), so
unlike the exact construction this one does respond to ties.

!!! note "The continuity correction is a choice"
    This package applies the ``1/2``, and so does R's `wilcox.test`, whose `correct`
    argument defaults to `TRUE`. The figures below say what dropping it would cost, not
    what either does.

    Without it the index is anticonservative. Across the 66 one-sample sizes ``n = 5:70``
    the interval comes out narrower than the exact one on 45 of them at
    ``1-\alpha = 0.90`` and on 7 at ``0.95``; with the correction, on 10 and on none.

The attained coverage of an approximate interval is not computed and is not guaranteed to
reach the nominal level.

### 6.5 One-sided intervals

A one-sided bound at level ``L`` is the corresponding endpoint of the two-sided interval
at level ``2L - 1``; that is, the two-sided ``\alpha`` used is ``2(1-L)`` rather than
``1-L``. The other endpoint is infinite.

Which of the two is kept follows the alternative the tail names, as in
[§2.1](@ref "2.1 Model, estimand, statistic"): the alternative ``\theta < 0`` is compatible
with an upper bound, so

```math
\text{left} \;\longrightarrow\; \bigl(-\infty,\, V_{(m-k)}\bigr) ,
\qquad
\text{right} \;\longrightarrow\; \bigl(V_{(k+1)},\, \infty\bigr) .
```

The endpoint kept is then the acceptance limit of the one-sided test of the same name, so
`pvalue` and `confint` given the same tail describe the same alternative. This is the
convention of [§7.1](@ref "7.1 One-sided intervals") of [The t-tests](@ref), of every other
test in this package that takes a tail, and of R under `alternative = "less"` and
`"greater"`. These four tests returned the other endpoint until
[issue #368](https://github.com/JuliaStats/HypothesisTests.jl/issues/368).

### 6.6 Zeros, ties, and degeneracy

**Zeros.** By [§2.1](@ref "2.1 Model, estimand, statistic") the one-sample statistic is
computed from the ``n`` non-zero differences. The pairwise estimates must be formed from the
same ``n`` observations, or the p-value and the interval describe different samples. For a
20-point sample containing five zeros there are 120 of them (``m = n(n+1)/2`` with
``n = 15``), not 210 (which is ``N(N+1)/2`` with ``N = 20``, the count obtained by retaining
the zeros). If every difference is zero, every pairwise estimate is zero and the interval
degenerates to the point ``0``.

**Ties on the exact interval.** [§6.3](@ref "6.3 Exact index") inverts the lattice
distribution of [§2.2.1](@ref "2.2.1 The lattice distribution") or
[§3.2.1](@ref "3.2.1 The lattice distribution"). Under ties the relevant null
distribution is the permutation one of
[§2.2.2](@ref "2.2.2 The permutation distribution") or
[§3.2.2](@ref "3.2.2 The permutation distribution"), so the attained coverage is approximate rather
than exact. With the endpoints included it errs above the nominal level; excluded, it can
err well below, through the endpoint atoms
[§6.2.1](@ref "6.2.1 What happens at the endpoints") measures.

The classical construction retains the untied distribution, and this page specifies it.
Two alternatives exist. One is to decline an exact interval under ties and fall back to
[§6.4](@ref "6.4 Normal-approximation index"), which is what R's `wilcox.test` does. The
other is to invert the tied permutation distribution itself, recomputing it at every
candidate shift, which is what `exactRankTests::wilcox.exact` does; being a different
construction it returns different endpoints, ``(-0.5, 1)`` at ``0.95`` on the sample of
[§8.2](@ref "8.2 One sample, five zeros and ties among the rest") where the classical
construction gives ``(-0.25, 0.75)``, and its reported estimate, ``0.375`` there, is not
the Hodges-Lehmann ``0.5`` either.

**Degeneracy.** If ``\mathbb{P}(W \le 0) > \alpha/2`` then no ``k`` of this form reaches the
nominal level, and the widest available interval, ``(V_{(1)}, V_{(m)})``, falls short of
it. Its coverage is ``1 - 2^{1-n}`` in the untied one-sample case: ``0.75`` at ``n = 3``,
``0.875`` at ``n = 4``, ``0.9375`` at ``n = 5``, so a ``0.95`` interval is unattainable
below ``n = 6`` and a ``0.99`` interval below ``n = 8``. At ``\mathbb{P}(W \le 0) = \alpha/2``
exactly there is nothing to warn about: ``k = 0`` attains the request exactly, as at
``n = 5`` and ``1 - \alpha = 0.9375``.

Returning a short interval as though it met the request would misstate the coverage, so
this package returns it and warns. On the exact route the warning names the coverage that
is attainable; for a one-sided bound it names the one-sided request and the one-sided
attainable coverage ``1 - \mathbb{P}(W \le 0)``, since only one tail can miss.

R is quieter here. It returns the same short interval with its `conf.level` attribute still
claiming the level that was asked for, and warns, renaming that attribute to the coverage
achieved, only once the shortfall exceeds ``\alpha/2``: silent at ``n = 5`` and
``1-\alpha = 0.95``, warning at two observations against two, where the attribute becomes
``2/3``.

The approximate route warns for a related but weaker reason. Its index rule
([§6.4](@ref "6.4 Normal-approximation index")) can ask for an order statistic outside the
pairwise estimates, which is not the same as the level being out of reach: at ``n = 8`` and
``1-\alpha = 0.99`` the rule asks for ``k = -1`` while the exact route attains ``0.9922``
and so meets the request. That warning therefore reports what happened, and says the exact
route is the way on, without asserting anything about coverage the route cannot compute.

Either way the interval returned is the widest the form admits, since it is the best there
is.

## 7. Relation to other implementations

R's `wilcox.test(conf.int = TRUE)` is the usual reference. Where it takes its exact route
it implements [§6.3](@ref "6.3 Exact index"), via the algorithm of
[Bauer (1972)](@cite bauer1972), and agrees digit for digit with this specification.

Three deliberate differences. The first two concern its approximate route:

  - It does not return order statistics. It solves numerically, with `uniroot`, for the
    shift at which the statistic crosses its critical value, so its endpoints lie near a
    pairwise estimate without being one: `3.0500354` where
    [§6.4](@ref "6.4 Normal-approximation index") gives `3.05`.
  - It continuity-corrects the interval, as
    [§6.4](@ref "6.4 Normal-approximation index") does, but not the point estimate, which
    it also root-finds rather than taking as the median of the pairwise estimates. Its reported
    estimate therefore drifts from ``\hat\theta``: `9.71184` against `9.675` on the sample
    of [§8.1](@ref "8.1 One sample, no ties and no zeros").
  - Under ties it declines an exact interval entirely and falls back to its approximate
    route, where [§6.6](@ref "6.6 Zeros, ties, and degeneracy") retains the classical
    construction instead.

It is also quieter in the degenerate case of
[§6.6](@ref "6.6 Zeros, ties, and degeneracy"), where it substitutes a lower-coverage
interval silently unless the shortfall exceeds ``\alpha/2``.

Its route-selection rule differs too. R takes the exact route when the sample is under 50,
each of the two for the rank sum test, and there are no ties, and for the signed rank test no
zeros either; this package instead lowers the size threshold under ties rather than abandoning
the exact route, as [§2.2](@ref "2.2 Null distribution of the signed rank statistic") and
[§3.2](@ref "3.2 Null distribution of the Mann-Whitney statistic") set out. Comparisons
against R must therefore set its `exact` argument explicitly, or the two implementations may
be running different constructions.

## 8. Worked values for the rank tests

Conformance vectors. Values are exact as printed unless a tolerance is implied by the
digits shown. The sessions are run when this page is built, so what is shown is what the
package returns; the tables carry the quantities behind them that no printed output shows.

Two labels in that printed output need decoding. `rank sums:` on the signed rank tests
gives ``W^+`` of [§2.1](@ref "2.1 Model, estimand, statistic") beside the midranks carried
by the negative observations, the two summing to ``n(n+1)/2``. And on the approximate
tests, `normal approximation (μ, σ):` reports the *centred* statistic of
[§2.2.3](@ref "2.2.3 The normal approximation") or
[§3.2.3](@ref "3.2.3 The normal approximation") beside the tie-corrected standard deviation,
so its first entry is ``W^+ - n(n+1)/4`` or ``U - n_x n_y/2``, not the null mean. The
other labels say what they mean: `Wilcoxon signed rank statistic:` is ``W^+`` itself, and
the Mann-Whitney tests report as `Location shift` the ``\Delta`` of
[§3.1](@ref "3.1 Model, estimand, statistic"), estimated by the Hodges-Lehmann median of
[§5](@ref "5. Point estimation").

### 8.1 One sample, no ties and no zeros

``n = 15``, ``d`` =
`[-7.8, -6.9, -4.7, 3.7, 6.5, 8.7, 9.1, 10.1, 10.8, 13.6, 14.4, 16.6, 20.2, 22.4, 23.5]`,
``m = 120``.

```jldoctest rank1
julia> using HypothesisTests

julia> d = [-7.8, -6.9, -4.7, 3.7, 6.5, 8.7, 9.1, 10.1, 10.8, 13.6, 14.4, 16.6, 20.2, 22.4, 23.5];

julia> t = SignedRankTest(d)
Exact Wilcoxon signed rank test
-------------------------------
Population details:
    parameter of interest:   Location parameter (pseudomedian)
    value under h_0:         0
    point estimate:          9.675
    95% confidence interval: (3.3, 15.5)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0034

Details:
    number of observations:         15
    non-zero observations:          15
    Wilcoxon signed rank statistic: 109.0
    rank sums:                      [109.0, 11.0]
    adjustment for ties:            0.0


julia> pvalue(t)
0.00335693359375

julia> hodgeslehmann(t)
9.675

julia> confint(t)
(3.3000000000000003, 15.5)

julia> confint(t; level = 0.90)
(4.449999999999999, 14.45)

julia> confint(t; level = 0.95, tail = :left)
(-Inf, 14.45)

julia> confint(t; level = 0.95, tail = :right)
(4.449999999999999, Inf)

julia> confint(ApproximateSignedRankTest(d))
(3.0500000000000003, 15.5)
```

The summary rounds the endpoints; `confint` returns them as computed, which is why the
tables below and the tests in this package compare to a tolerance.

| quantity | value |
|:---|:---|
| ``\hat\theta`` ([§5](@ref "5. Point estimation")) | `9.675` |
| median of ``d``, for comparison | `10.1` |
| exact index ([§6.3](@ref "6.3 Exact index")) at ``1-\alpha = 0.95`` | ``k = 25``, ``C_\alpha = 26`` |
| attained coverage | `0.95209`; at ``k=26`` it is `0.94464` |
| approximate index ([§6.4](@ref "6.4 Normal-approximation index")) at ``1-\alpha = 0.95`` | ``\sigma = 17.60682``, ``C_\alpha = 25`` |

The one-sided calls illustrate [§6.5](@ref "6.5 One-sided intervals"): each bound at
``0.95`` is an endpoint of the two-sided interval at ``0.90``, the upper one for
`tail = :left` and the lower for `tail = :right`. The last call shows the approximate route
moving the lower endpoint one order statistic outwards, from ``V_{(26)} = 3.3`` to
``V_{(25)} = 3.05``; its upper endpoint is ``15.5`` either way, two of the Walsh averages
there being equal.

### 8.2 One sample, five zeros and ties among the rest

``N = 20``, ``d`` = `[0, 0, 0, 0.5, 0.5, 1, -0.5, -1, 1.5, -1.5, 0.5, 0, 1, -0.5, 2, 0, 0.5, -1, 1, 0.5]`,
so ``n = 15`` and ``T(|d|) = 462``. Both counts are reported, as
[§2.1](@ref "2.1 Model, estimand, statistic") asks, and it is ``n = 15`` that everything is
computed from.

```jldoctest rank2
julia> using HypothesisTests

julia> d = [0, 0, 0, 0.5, 0.5, 1, -0.5, -1, 1.5, -1.5, 0.5, 0, 1, -0.5, 2, 0, 0.5, -1, 1, 0.5];

julia> t = SignedRankTest(d)
Exact Wilcoxon signed rank test
-------------------------------
Population details:
    parameter of interest:   Location parameter (pseudomedian)
    value under h_0:         0
    point estimate:          0.5
    95% confidence interval: (-0.25, 0.75)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.3072

Details:
    number of observations:         20
    non-zero observations:          15
    Wilcoxon signed rank statistic: 78.5
    rank sums:                      [78.5, 41.5]
    adjustment for ties:            462.0


julia> pvalue(t)
0.30718994140625

julia> confint(t; level = 0.90)
(-0.25, 0.75)

julia> length(HypothesisTests.walsh_averages(filter(!iszero, d)))
120
```

| quantity | value |
|:---|:---|
| ``\hat\theta`` | `0.5` |
| number of pairwise estimates ([§6.6](@ref "6.6 Zeros, ties, and degeneracy")) | `120` = ``15 \cdot 16 / 2``, against `210` = ``20 \cdot 21 / 2`` if zeros were retained |
| two-sided p-value ([§2.3](@ref "2.3 p-values"), tied branch) | `0.30719`, by the enumeration of [§2.2.2](@ref "2.2.2 The permutation distribution") |

### 8.3 Two samples, no ties

``x`` = `1:10`, ``y`` = `2.1, 4.1, …, 20.1`; ``m = 100``, ``U = 20``.

```jldoctest rank3
julia> using HypothesisTests

julia> x = collect(1.0:10.0); y = collect(2.1:2.0:20.1);

julia> t = MannWhitneyUTest(x, y)
Exact Mann-Whitney U test
-------------------------
Population details:
    parameter of interest:   Location shift
    value under h_0:         0
    point estimate:          -5.6
    95% confidence interval: (-11.1, -0.1)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0232

Details:
    number of observations in each group: [10, 10]
    Mann-Whitney-U statistic:             20.0
    rank sums:                            [75.0, 135.0]
    adjustment for ties:                  0.0


julia> pvalue(t)
0.02323063932971052

julia> hodgeslehmann(t)
-5.6

julia> confint(t)
(-11.1, -0.09999999999999964)
```

| quantity | value |
|:---|:---|
| ``\hat\Delta`` | `-5.6` |
| exact index at ``1-\alpha = 0.95`` | ``k = 23``, ``C_\alpha = 24`` |
| ``\mathbb{P}(U \le 23)``, ``\mathbb{P}(U \le 24)`` | `0.021629`, `0.026213` |
| attained coverage | `0.95674` |

### 8.4 Two samples, ties

``x`` = `1:10`, ``y`` = `2, 4, …, 24`, so ``N = 22`` and five values are tied across the
samples. Here ``N > 10``, so the routing of
[§3.2](@ref "3.2 Null distribution of the Mann-Whitney statistic") sends tied data this size
to the approximate test, and `method = :exact` is what reaches the exact one.

```jldoctest rank4
julia> using HypothesisTests

julia> x = collect(1.0:10.0); y = collect(2.0:2.0:24.0);

julia> t = MannWhitneyUTest(x, y)
Approximate Mann-Whitney U test
-------------------------------
Population details:
    parameter of interest:   Location shift
    value under h_0:         0
    point estimate:          -7.5
    95% confidence interval: (-14.0, -1.0)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0146

Details:
    number of observations in each group: [10, 12]
    Mann-Whitney-U statistic:             22.5
    rank sums:                            [77.5, 175.5]
    adjustment for ties:                  30.0
    normal approximation (μ, σ):          (-37.5, 15.1443)


julia> confint(MannWhitneyUTest(x, y; method = :exact))
(-14.0, -1.0)

julia> confint(MannWhitneyUTest(y, x; method = :exact))
(1.0, 14.0)
```

On this sample the two routes agree: the normal index lands on the same pair of order
statistics. R declines an exact interval under ties altogether
([§7](@ref "7. Relation to other implementations")) and solves numerically on its approximate
route, which lands here too, reporting `(-13.99997, -1.00003)` for the same pair. The last
call is the sample exchange of [§3.1](@ref "3.1 Model, estimand, statistic"): the interval
negates and reverses, and the two-sided p-value does not move.

## 9. The rank tests in this package

Four types implement the two procedures: [`ExactSignedRankTest`](@ref) and
[`ApproximateSignedRankTest`](@ref) for one sample,
[`ExactMannWhitneyUTest`](@ref) and [`ApproximateMannWhitneyUTest`](@ref) for two. The
`Exact*` types implement [§6.3](@ref "6.3 Exact index") and the `Approximate*` types
[§6.4](@ref "6.4 Normal-approximation index").

[`SignedRankTest`](@ref) and [`MannWhitneyUTest`](@ref) are functions rather than types.
Each ranks the sample, applies the selection rule of
[§2.2](@ref "2.2 Null distribution of the signed rank statistic") or
[§3.2](@ref "3.2 Null distribution of the Mann-Whitney statistic") unless the `method`
keyword says otherwise, and returns one of the four types above. There is no type of either
name to dispatch on or annotate with, and no supertype joining a procedure's two.

`pvalue` implements [§2.3](@ref "2.3 p-values") and [§3.3](@ref "3.3 p-values"), `confint`
implements [§6](@ref "6. Interval estimation"), and [`hodgeslehmann`](@ref) implements
[§5](@ref "5. Point estimation").

The exact null distributions of [§2.2.1](@ref "2.2.1 The lattice distribution") and
[§3.2.1](@ref "3.2.1 The lattice distribution") come from StatsFuns, which carries the
normalisation inside its recursions rather than accumulating lattice counts, so it takes the
first of the two routes those sections describe. The `[compat]` floor of StatsFuns 2.2.1 is
what makes the numerical care of both a settled question here rather than a caveat. The tied
routes are this package's own, and are bounded rather than corrected.

**Departures from this specification**, recorded here because
[§7](@ref "7. Relation to other implementations") asks the same of other implementations.

  - Every test carries a `median` field: the sample median for the signed rank tests, the
    difference of sample medians for the Mann-Whitney ones. Both are descriptive statistics
    rather than the estimand of [§5](@ref "5. Point estimation"), which is reported separately
    and is what the interval is built around. The signed rank one is taken over all ``N``
    differences, zeros included, so it is the one number these types carry that the discard
    of [§2.1](@ref "2.1 Model, estimand, statistic") does not reach. It is also that field,
    rather than anything the procedure needs, that makes an empty sample throw instead of
    returning the degenerate p-value of [§3.3](@ref "3.3 p-values").
  - The pairwise estimates are materialised as written in
    [§4.1](@ref "4.1 Definitions") rather than found by selection as that section describes,
    which bounds the usable sample size at
    [`MAX_PAIRWISE_ESTIMATES`](@ref HypothesisTests.MAX_PAIRWISE_ESTIMATES) of them.
  - The tied enumerations of [§2.2.2](@ref "2.2.2 The permutation distribution") and
    [§3.2.2](@ref "3.2.2 The permutation distribution") are bounded by
    [`MAX_EXACT_ENUMERATION_N`](@ref HypothesisTests.MAX_EXACT_ENUMERATION_N), past which the
    p-value is refused rather than computed.
  - Inverting the exact null distribution of [§6.3](@ref "6.3 Exact index") runs a lattice
    recursion per bisection step, which the definition there does not charge for. The
    two-sample recursion is the expensive one and is bounded separately, below the
    materialisation above, by
    [`MAX_EXACT_CI_ESTIMATES`](@ref HypothesisTests.MAX_EXACT_CI_ESTIMATES): past that the
    exact interval is refused, and the normal-approximation interval of
    [§6.4](@ref "6.4 Normal-approximation index"), which inverts a closed form, is not. The
    one-sample recursion is cheap enough that the materialisation bound is the only one its
    exact interval meets.
  - Under ties an `Exact*` test enumerates for its p-value but inverts the untied lattice for
    its interval, so the two disagree about ties. That is the classical construction, and
    [§6.6](@ref "6.6 Zeros, ties, and degeneracy") says which way it errs.

## 10. References

Both procedures originate with [Wilcoxon (1945)](@cite wilcoxon1945), which introduced the
signed rank and the rank sum tests together; [Mann and Whitney (1947)](@cite mann1947)
arrived at the two-sample test independently and in the counting form, which
[§3](@ref "3. The two-sample procedure (Wilcoxon rank sum, Mann-Whitney U)") reconciles
with the rank sum.
The estimator is due to [Hodges and Lehmann (1963)](@cite hodges1963). The two-sample
interval of [§6](@ref "6. Interval estimation") is commonly called the **Moses interval**,
after the chapter L. E. Moses contributed to Walker and Lev's *Statistical Inference*
(1953); the one-sample counterpart is usually credited to Tukey. Both constructions, with
tables, are in [Hollander and Wolfe (1973)](@cite hollander1973) at pages 27–33 and 68–75.
[Bauer (1972)](@cite bauer1972) gives the order-statistic algorithm.
