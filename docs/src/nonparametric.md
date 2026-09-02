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

Both are built on a set of pairwise estimates: the Walsh averages ``(dᵢ + dⱼ)/2`` for the
signed rank tests, the cross-group differences ``xᵢ - yⱼ`` for the Mann-Whitney tests.
`confint` returns a pair of order statistics of that set, chosen by inverting the test,
and `hodgeslehmann` returns its median.

```@docs
hodgeslehmann
```

The Hodges-Lehmann estimate is the point estimate the interval brackets. It is generally
not the sample median: exact symmetry makes the two agree, but they can also agree by
coincidence, as on `[0, 2, 2, 7]`, where both are `2`. Agreement is therefore no evidence
of symmetry.

Three named bounds apply, and past any of them the tests raise rather than run unbounded.
The first bounds the tied-data enumeration a p-value runs, and `method = :approximate` is
the way past it. The second bounds the set of pairwise estimates, which `confint` and
`hodgeslehmann` form whichever route the p-value took, so `method` is no help there; it is a
bound on memory, which is all the approximate interval spends. The third bounds the exact
interval, which spends time as well, a lattice recursion for every candidate endpoint it
considers, and so is bounded far below the second.

All three raise `ComputationTooLarge`, which is its own type rather than an `ArgumentError`
so that `show` can drop a refused interval line without also hiding real errors.

```@docs
HypothesisTests.MAX_EXACT_ENUMERATION_N
HypothesisTests.MAX_PAIRWISE_ESTIMATES
HypothesisTests.MAX_EXACT_CI_ESTIMATES
HypothesisTests.ComputationTooLarge
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
