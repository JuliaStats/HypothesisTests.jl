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

__precompile__()

module HypothesisTests

using Distributions, Roots, StatsBase, Compat
using Combinatorics: combinations
using Rmath: pwilcox, psignrank

import StatsBase.confint

export testname, pvalue, confint
@compat abstract type HypothesisTest end

check_same_length(x::AbstractVector, y::AbstractVector) = if length(x) != length(y)
    throw(DimensionMismatch("Vectors must be the same length"))
end

# Basic function for finding a p-value given a distribution and tail
pvalue(dist::ContinuousUnivariateDistribution, x::Number; tail=:both) =
    if tail == :both
        min(2 * min(cdf(dist, x), ccdf(dist, x)), 1.0)
    elseif tail == :left
        cdf(dist, x)
    elseif tail == :right
        ccdf(dist, x)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end

pvalue(dist::DiscreteUnivariateDistribution, x::Number; tail=:both) =
    if tail == :both
        min(2 * min(ccdf(dist, x-1), cdf(dist, x)), 1.0)
    elseif tail == :left
        cdf(dist, x)
    elseif tail == :right
        ccdf(dist, x-1)
    else
        throw(ArgumentError("tail=$(tail) is invalid"))
    end

function check_alpha(alpha::Float64)
    if alpha <= 0 || alpha >= 0.5
        throw(ArgumentError("alpha $alpha not in range (0, 0.5)"))
    end
end

# Pretty-print
function Base.show{T<:HypothesisTest}(io::IO, test::T)
    println(io, testname(test))
    println(io, repeat("-", length(testname(test))))

    # population details
    has_ci = applicable(StatsBase.confint, test)
    (param_name, param_under_h0, param_estimate) = population_param_of_interest(test)
    println(io, "Population details:")
    println(io, "    parameter of interest:   $param_name")
    println(io, "    value under h_0:         $param_under_h0")
    println(io, "    point estimate:          $param_estimate")
    if has_ci
        println(io, "    95% confidence interval: $(StatsBase.confint(test))")
    end
    println(io)

    # test summary
    p = pvalue(test)
    outcome = if p > 0.05 "fail to reject" else "reject" end
    tail = default_tail(test)
    println(io, "Test summary:")
    println(io, "    outcome with 95% confidence: $outcome h_0")
    if tail == :both
        println(io, "    two-sided p-value:           $p")
    elseif tail == :left || tail == :right
        println(io, "    one-sided p-value:           $p")
    else
        println(io, "    p-value:                     $p")
    end
    println(io)

    # further details
    println(io, "Details:")
    show_params(io, test, "    ")
end

# parameter of interest: name, value under h0, point estimate
population_param_of_interest{T<:HypothesisTest}(test::T) = ("not implemented yet", NaN, NaN)

# is the test one- or two-sided
default_tail(test::HypothesisTest) = :undefined

function show_params{T<:HypothesisTest}(io::IO, test::T, ident="")
    fieldidx = find(Bool[t<:Number for t in T.types])
    if !isempty(fieldidx)
        lengths = [length(string(T.names[i])) for i in fieldidx]
        maxlen = maximum(lengths)

        for i = 1:length(fieldidx)
            name = T.names[fieldidx[i]]
            println(io, ident, repeat(" ", maxlen-lengths[i]),
                      replace(string(name), "_", " "),
                      " = $(getfield(test, name))")
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
end
