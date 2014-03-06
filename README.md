HypothesisTests.jl
===========

[![Build Status](https://travis-ci.org/JuliaStats/HypothesisTests.jl.png?branch=master)](https://travis-ci.org/JuliaStats/HypothesisTests.jl)

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
```

## Documentation

Full documentation available at [Read the Docs](http://hypothesistestsjl.readthedocs.org/en/latest/index.html).

## Credits

Original API suggested by [John Myles White](https://github.com/johnmyleswhite).
