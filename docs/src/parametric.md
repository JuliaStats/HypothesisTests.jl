# Parametric tests

```@setup ht
using HypothesisTests
```

## Power divergence test
```@docs
PowerDivergenceTest
```

## Pearson chi-squared test
```@docs
ChisqTest
```

## Multinomial likelihood ratio test
```@docs
MultinomialLRTest
```

## t-test
```@docs
OneSampleTTest
EqualVarianceTTest
UnequalVarianceTTest
```

## z-test
```@docs
OneSampleZTest
EqualVarianceZTest
UnequalVarianceZTest
```

## F-test
```@docs
VarianceFTest
```

## One-way ANOVA Test

```@docs
OneWayANOVATest
```

Example:

```@example ht
groups = [
    [6, 8, 4, 5, 3, 4],
    [8, 12, 9, 11, 6, 8],
    [13, 9, 11, 8, 7, 12]
]
t = OneWayANOVATest(groups...)
show(IOContext(stdout, :table => true), t)
```

## Levene's Test

```@docs
LeveneTest
```

## Brown-Forsythe Test

```@docs
BrownForsytheTest
```
