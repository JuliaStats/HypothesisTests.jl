# Nonparametric tests

## Anderson-Darling test

Available are both one-sample and ``k``-sample tests.

```@docs
OneSampleADTest
KSampleADTest
```

## Binomial test

```@docs
BinomialTest
confint(::BinomialTest)
```

## Fisher exact test

```@docs
FisherExactTest
confint(::FisherExactTest)
pvalue(::FisherExactTest)
```

## Kolmogorov-Smirnov test

Available are an exact one-sample test and approximate (i.e. asymptotic) one- and two-sample
tests.

```@docs
ExactOneSampleKSTest
ApproximateOneSampleKSTest
ApproximateTwoSampleKSTest
```

## Kruskal-Wallis rank sum test

```@docs
KruskalWallisTest
```

## Mann-Whitney U test

```@docs
MannWhitneyUTest
ExactMannWhitneyUTest
ApproximateMannWhitneyUTest
```

## Sign test

```@docs
SignTest
```

## Wald-Wolfowitz independence test

```@docs
WaldWolfowitzTest
```

## Wilcoxon signed rank test

```@docs
SignedRankTest
ExactSignedRankTest
ApproximateSignedRankTest
```

## Rank tests: intervals and point estimates

The Wilcoxon signed rank and Mann-Whitney U tests share their interval and estimator
machinery.

`confint` inverts the test to a pair of order statistics of a contrast set — the Walsh
averages ``(dᵢ + dⱼ)/2`` for the signed rank tests, the cross-group differences
``xᵢ - yⱼ`` for the Mann-Whitney tests. Which order statistics depends on the test: the
`Exact*` types invert the exact null distribution, conservatively, and the
`Approximate*` types invert the normal approximation with a continuity correction and,
where relevant, a tie-corrected variance. Both reproduce R's `wilcox.test` for the
corresponding setting of its `exact` argument.

`hodgeslehmann` returns the median of that same contrast set. This is the point estimate
the interval brackets, and the one these tests report as their parameter of interest; it
is generally not the sample median.

```@docs
hodgeslehmann
```

Two bounds apply to the exact route, both because the underlying `StatsFuns` routines
accumulate lattice counts in `Int` and go silently wrong once those overflow. Past either
bound the exact tests raise rather than return an invalid answer, and `method =
:approximate` is the way on.

```@docs
HypothesisTests.MAX_EXACT_SIGNRANK_N
HypothesisTests.MAX_EXACT_WILCOX_BINOMIAL
HypothesisTests.MAX_EXACT_ENUMERATION_N
```

## Permutation test

```@docs
ExactPermutationTest
ApproximatePermutationTest
```

## Fligner-Killeen test

```@docs
FlignerKilleenTest
```

## Shapiro-Wilk test
```@docs
ShapiroWilkTest
```
