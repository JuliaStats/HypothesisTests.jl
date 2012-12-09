Wilcoxon.jl
===========

### Wilcoxon signed rank and rank sum hypothesis tests

This package implements the Wilcoxon signed rank test and rank sum (Mann-Whitney U) test in Julia.

The rank sum test (```mannwhitneyu```) uses the ```pwilcox``` function from [Rmath](http://www.r-bloggers.com/julia-functions-for-the-rmath-library/) to compute exact p-values when there are <= 50 samples and no tied ranks. If there are tied ranks and <= 10 samples, the test computes the exact p value by enumerating all possible combinations. In all other cases, the test generates a p value using the normal approximation, adjusting the variance for ties.

The signed rank test (```signrank```) uses the psignrank function from [Rmath](http://www.r-bloggers.com/julia-functions-for-the-rmath-library/) to compute exact p-values when there are <= 50 samples and no tied ranks. For larger samples, or samples with tied ranks, the function uses the normal approximation, adjusting the variance for ties.