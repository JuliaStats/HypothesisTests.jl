HypothesisTests.jl
===========

[![Build Status](https://travis-ci.org/JuliaStats/HypothesisTests.jl.svg?branch=master)](https://travis-ci.org/JuliaStats/HypothesisTests.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaStats/HypothesisTests.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaStats/HypothesisTests.jl?branch=master)
[![HypothesisTests](http://pkg.julialang.org/badges/HypothesisTests_0.5.svg)](http://pkg.julialang.org/?pkg=HypothesisTests)
[![HypothesisTests](http://pkg.julialang.org/badges/HypothesisTests_0.6.svg)](http://pkg.julialang.org/?pkg=HypothesisTests)

This package implements several hypothesis tests in Julia.

## Quick start

Some examples:

```julia
using HypothesisTests

pvalue(OneSampleTTest(x))
pvalue(OneSampleTTest(x), tail=:left)
pvalue(OneSampleTTest(x), tail=:right)
confint(OneSampleTTest(x))
confint(OneSampleTTest(x, tail=:left))
confint(OneSampleTTest(x, tail=:right))
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

Full documentation available at
[https://JuliaStats.github.io/HypothesisTests.jl/latest/
](https://JuliaStats.github.io/HypothesisTests.jl/latest/).
