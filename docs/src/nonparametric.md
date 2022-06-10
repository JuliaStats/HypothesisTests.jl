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

## Permutation test

```@docs
ExactPermutationTest
ApproximatePermutationTest
```

## Fligner-Killeen test

```@docs
FlignerKilleenTest
```
