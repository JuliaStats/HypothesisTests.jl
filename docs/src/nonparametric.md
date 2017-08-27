# Nonparametric tests

## Anderson-Darling test

The null hypothesis of the Anderson–Darling test is that a dataset comes from a certain
distribution; the reference distribution can be specified explicitly (one-sample test).
``k``-sample Anderson–Darling tests are available for testing whether several samples are
coming from a single population drawn from the distribution function which does not have to
be specified.

```@docs
OneSampleADTest
KSampleADTest
```

## Binomial test

```@docs
BinomialTest
```

## Fisher exact test

```@docs
FisherExactTest
```

## Kolmogorov–Smirnov test

The null hypothesis of the Kolmogorov–Smirnov test is that a dataset comes from a certain
distribution; the reference distribution can be specified explicitly (one-sample test) or by
an empirical sample (two-sample test). The alternative hypothesis is that the cumulative
distributions of the sample is different (`tail = :both`: default), smaller (`tail = :left`),
or larger (`tail = :right`) than the reference cumulative distribution. The exact test is
based on the exact distribution of the differences whereas the approximate test is derived
from its asymptotic distribution.

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

## Wilcoxon signed rank test

```@docs
SignedRankTest
ExactSignedRankTest
ApproximateSignedRankTest
```
