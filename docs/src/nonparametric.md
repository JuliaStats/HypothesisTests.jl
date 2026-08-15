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

Both are built on a set of pairwise contrasts: the Walsh averages ``(dᵢ + dⱼ)/2`` for the
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

Two named bounds apply, one to the tied-data enumeration and one to the contrast set.
Past either the tests raise rather than run unbounded, and `method = :approximate` is the
way on.

```@docs
HypothesisTests.MAX_EXACT_ENUMERATION_N
HypothesisTests.MAX_RANK_CONTRASTS
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
