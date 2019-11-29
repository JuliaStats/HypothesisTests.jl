## HypothesisTests.jl

*HypothesisTests.jl* is a Julia package that implements a wide range of hypothesis tests.

- **Build & Testing Status:**
  [![Build Status](https://travis-ci.org/JuliaStats/HypothesisTests.jl.svg?branch=master)
  ](https://travis-ci.org/JuliaStats/HypothesisTests.jl)
  [![Coverage Status](https://coveralls.io/repos/JuliaStats/HypothesisTests.jl/badge.svg?branch=master)
  ](https://coveralls.io/r/JuliaStats/HypothesisTests.jl?branch=master)

- **Documentation**: [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://JuliaStats.github.io/HypothesisTests.jl/latest/

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
