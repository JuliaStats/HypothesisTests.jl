HypothesisTests.jl
===========

This package implements several hypothesis tests in Julia.

## Quick start

To compute p-values:

```julia
p_value(MannWhitneyUTest, x, y)
p_value(SignedRankTest, x, y)
p_value(SignedRankTest, x)
p_value(RayleighTest, complex_numbers)
p_value(RayleighTest, angles_in_radians)
left_p_value(MannWhitneyUTest, x, y)
right_p_value(SignedRankTest, x, y)
```

To compute test statistics:

```julia
test_statistic(MannWhitneyUTest, x, y)
```

To show test statistics and p-values:

```julia
MannWhitneyUTest(x, y)
SignedRankTest(x, y)
RayleighTest(x)
```

## Implementation notes

```MannWhitneyUTest``` uses the ```pwilcox``` function from [Rmath](http://www.r-bloggers.com/julia-functions-for-the-rmath-library/) to compute exact p-values when there are <= 50 samples and no tied ranks. If there are tied ranks and <= 10 samples, the test computes the exact p value by enumerating all possible combinations. In all other cases, the test generates a p value using the normal approximation, adjusting the variance for ties. Behavior can be controlled by using ```ExactMannWhitneyUTest``` and ```ApproximateMannWhitneyUTest```, although the former is very slow for even moderately sized samples with tied ranks.

```SignedRankTest``` uses the ```psignrank``` function from [Rmath](http://www.r-bloggers.com/julia-functions-for-the-rmath-library/) to compute exact p-values when there are <= 50 samples and no tied ranks. If there are tied ranks and <= 15 samples, the test computes the exact p value by enumerating all possible combinations of signs. For larger samples, or samples with tied ranks, the function uses the normal approximation, adjusting the variance for ties. Behavior can be controlled by using ```ExactSignedRankTest``` and ```ApproximateSignedRankTest```, although the former is very slow for even moderately sized samples with tied ranks.

```RayleighTest``` uses the formula from Fisher, 1993 to compute the p-value.

## Credits

API suggested by [John Myles White](https://github.com/johnmyleswhite).