## HypothesisTests.jl

*HypothesisTests.jl* is a Julia package that implements a wide range of hypothesis tests.

- **Build & Testing Status:**
[![Build Status](https://github.com/JuliaStats/HypothesisTests.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/JuliaStats/HypothesisTests.jl/actions/workflows/ci.yml?query=branch%3Amaster)
[![Coverage Status](https://codecov.io/gh/JuliaStats/HypothesisTests.jl/branch/master/graph/badge.svg?token=ztSoTXYVhb)](https://codecov.io/gh/JuliaStats/HypothesisTests.jl)

- **Documentation**: [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-latest-url]: http://JuliaStats.github.io/HypothesisTests.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://JuliaStats.github.io/HypothesisTests.jl/stable/

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
