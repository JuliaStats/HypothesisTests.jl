# HypothesisTests.jl
# Hypothesis tests in Julia
#
# Copyright (C) 2012   Simon Kornblith
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

module HypothesisTests

using Statistics, Random, LinearAlgebra
using Distributions, Roots, StatsBase
using Combinatorics: combinations, permutations
using Rmath: pwilcox, psignrank

import StatsBase.confint

export testname, pvalue, confint
abstract type HypothesisTest end

check_same_length(x::AbstractVector, y::AbstractVector) = if length(x) != length(y)
    throw(DimensionMismatch("Vectors must be the same length"))
end

"""
    confint(test::HypothesisTest; level = 0.95, tail = :both)

Compute a confidence interval C with coverage `level`.

If `tail` is `:both` (default), then a two-sided confidence interval is returned. If `tail`
is `:left` or `:right`, then a one-sided confidence interval is returned.

!!! note
    Most of the implemented confidence intervals are *strongly consistent*, that is, the
    confidence interval with coverage `level` does not contain the test statistic under
    ``h_0`` if and only if the corresponding test rejects the null hypothesis
    ``h_0: θ = θ_0``:
    ```math
        C (x, level) = \\{θ : p_θ (x) > 1 - level\\},
    ```
    where ``p_θ`` is the [`pvalue`](@ref) of the corresponding test.
"""
function confint end

"""
    pvalue(test::HypothesisTest; tail = :both)

Compute the p-value for a given significance test.

If `tail` is `:both` (default), then the p-value for the two-sided test is returned. If
`tail` is `:left` or `:right`, then a one-sided test is performed.
"""
function pvalue end

# Basic function for finding a p-value given a distribution and tail
function pvalue(dist::ContinuousUnivariateDistribution, x::Number; tail=:both)
    check_tail(tail)

    if tail == :both
        p = 2 * min(cdf(dist, x), ccdf(dist, x))
        min(p, oneunit(p)) # if P(X = x) > 0, then possibly p > 1
    elseif tail == :left
        cdf(dist, x)
    else # tail == :right
        ccdf(dist, x)
    end
end

function pvalue(dist::DiscreteUnivariateDistribution, x::Number; tail=:both)
    check_tail(tail)

    if tail == :both
        p = 2 * min(ccdf(dist, x-1), cdf(dist, x))
        min(p, oneunit(p)) # if P(X = x) > 0, then possibly p > 1
    elseif tail == :left
        cdf(dist, x)
    else # tail == :right
        ccdf(dist, x-1)
    end
end

function check_level(level::Float64)
    if level >= 1 || level <= 0.5
        throw(ArgumentError("coverage level $level not in range (0.5, 1)"))
    end
end

function check_tail(tail::Symbol)
    if tail !== :both && tail !== :left && tail !== :right
        throw(ArgumentError("tail=$(tail) is invalid"))
    end
end

# Pretty-print
function Base.show(_io::IO, test::T) where T<:HypothesisTest
    io = IOContext(_io, :compact=>get(_io, :compact, true))
    println(io, testname(test))
    println(io, repeat("-", length(testname(test))))

    # population details
    has_ci = applicable(StatsBase.confint, test)
    (param_name, param_under_h0, param_estimate) = population_param_of_interest(test)
    println(io, "Population details:")
    println(io, "    parameter of interest:   $param_name")
    print(io, "    value under h_0:         ")
    show(io, param_under_h0)
    println(io)
    print(io, "    point estimate:          ")
    show(io, param_estimate)
    println(io)

    if has_ci
        ci = map(x -> round.(x; sigdigits=4, base=10), StatsBase.confint(test))
        print(io, "    95% confidence interval: ")
        show(io, ci)
        println(io)
    end
    println(io)

    # test summary
    p = pvalue(test)
    outcome = if p > 0.05 "fail to reject" else "reject" end
    tail = default_tail(test)
    pval = StatsBase.PValue(p)
    println(io, "Test summary:")
    println(io, "    outcome with 95% confidence: $outcome h_0")
    if tail == :both
        println(io, "    two-sided p-value:           $pval")
    elseif tail == :left || tail == :right
        println(io, "    one-sided p-value:           $pval")
    else
        println(io, "    p-value:                     $pval")
    end
    println(io)

    # further details
    println(io, "Details:")
    show_params(io, test, "    ")
end

# parameter of interest: name, value under h0, point estimate
population_param_of_interest(test::T) where {T<:HypothesisTest} = ("not implemented yet", NaN, NaN)

# is the test one- or two-sided
default_tail(test::HypothesisTest) = :undefined

function show_params(io::IO, test::T, ident="") where T<:HypothesisTest
    fieldidx = findall(Bool[t<:Number for t in T.types])
    if !isempty(fieldidx)
        lengths = [length(string(T.names[i])) for i in fieldidx]
        maxlen = maximum(lengths)

        for i = 1:length(fieldidx)
            name = T.names[fieldidx[i]]
            print(io, ident, repeat(" ", maxlen-lengths[i]),
                      replace(string(name), "_", " ", " = "))
            show(io, getfield(test, name))
            println(io)
        end
    end
end

include("deprecated.jl")
include("common.jl")

include("binomial.jl")
include("circular.jl")
include("fisher.jl")
include("kolmogorov_smirnov.jl")
include("kruskal_wallis.jl")
include("mann_whitney.jl")
include("t.jl")
include("z.jl")
include("wilcoxon.jl")
include("power_divergence.jl")
include("anderson_darling.jl")
include("box_test.jl")
include("breusch_godfrey.jl")
include("augmented_dickey_fuller.jl")
include("jarque_bera.jl")
include("durbin_watson.jl")
include("permutation.jl")
include("hotelling.jl")
include("bartlett.jl")
include("wald_wolfowitz.jl")
include("f.jl")
include("correlation.jl")
include("diebold_mariano.jl")
include("clark_west.jl")
include("white.jl")
include("var_equality.jl")

end
