HypothesisTests.jl
===========

[![Build Status](https://travis-ci.org/simonster/HypothesisTests.jl.png?branch=master)](https://travis-ci.org/simonster/HypothesisTests.jl)

This package implements several hypothesis tests in Julia.

## Quick start

Some examples:

```julia
using HypothesisTests

pvalue(OneSampleTTest(x))
pvalue(OneSampleTTest(x), tail=:left)
pvalue(OneSampleTTest(x), tail=:right)
ci(OneSampleTTest(x))
ci(OneSampleTTest(x, tail=:left))
ci(OneSampleTTest(x, tail=:right))
OneSampleTTest(x).t
OneSampleTTest(x).df

pvalue(OneSampleTTest(x, y))
pvalue(EqualVarianceTTest(x, y))
pvalue(UnequalVarianceTTest(x, y))

pvalue(MannWhitneyUTest(x, y))
pvalue(SignedRankTest(x, y))
pvalue(SignedRankTest(x))

pvalue(RayleighTest(complex_numbers))
pvalue(RayleighTest(angles_in_radians))
```

## Implementation notes

```MannWhitneyUTest``` uses the ```pwilcox``` function from [Rmath](http://www.r-bloggers.com/julia-functions-for-the-rmath-library/) to compute exact p-values when there are <= 50 samples and no tied ranks. If there are tied ranks and <= 10 samples, the test computes the exact p value by enumerating all possible combinations. In all other cases, the test generates a p value using the normal approximation, adjusting the variance for ties. Behavior can be controlled by using ```ExactMannWhitneyUTest``` and ```ApproximateMannWhitneyUTest```, although the former is very slow for even moderately sized samples with tied ranks.

```SignedRankTest``` uses the ```psignrank``` function from [Rmath](http://www.r-bloggers.com/julia-functions-for-the-rmath-library/) to compute exact p-values when there are <= 50 samples and no tied ranks. If there are tied ranks and <= 15 samples, the test computes the exact p value by enumerating all possible combinations of signs. For larger samples, or samples with tied ranks, the function uses the normal approximation, adjusting the variance for ties. Behavior can be controlled by using ```ExactSignedRankTest``` and ```ApproximateSignedRankTest```, although the former is very slow for even moderately sized samples with tied ranks.

```RayleighTest``` uses the formula from Fisher, 1993 to compute the p-value.

## Credits

Original API suggested by [John Myles White](https://github.com/johnmyleswhite).