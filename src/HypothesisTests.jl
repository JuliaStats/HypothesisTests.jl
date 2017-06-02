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

# Utility to get an optional field (e.g., :tail and :alpha)
getfield(value::HypothesisTest, name::Symbol, default::Any) =
    if name in fieldnames(value)
        Base.getfield(value, name)
    else
        default
    end

get_tail(test::HypothesisTest) = getfield(test, :tail, default_tail(test))
get_alpha(test::HypothesisTest) = getfield(test, :alpha, 0.05)

# Utility for pretty-printing: Append white space so that length(with_trailing_whitespace(s)) = max(len, length(s))
with_trailing_whitespace(s::String, len::Int) = s * join(repeat([" "], outer=max(len - length(s), 0)), "")

# Pretty-print
function Base.show{T<:HypothesisTest}(io::IO, test::T)
    println(io, testname(test))
    println(io, repeat("-", length(testname(test))))
    
    conf_string = string((1 - get_alpha(test)) * 100) * "%"

    # population details
    has_ci = applicable(StatsBase.confint, test)
    label_len = has_ci ? length(conf_string) + 21 : 22 # 21 is length of ' confidence interval:', 22 that of 'parameter of interest:'
    (param_name, param_under_h0, param_estimate) = population_param_of_interest(test)
    println(io, "Population details:")
    println(io, "    $(with_trailing_whitespace("parameter of interest:", label_len)) $param_name")
    println(io, "    $(with_trailing_whitespace("value under h_0", label_len)) $param_under_h0")
    println(io, "    $(with_trailing_whitespace("point estimate", label_len)) $param_estimate")
    if has_ci
        println(io, "    $conf_string confidence interval: $(StatsBase.confint(test))")
    end
    println(io)

    # test summary
    tail = get_tail(test)
    p = pvalue(test, tail=tail)
    outcome = if p > get_alpha(test) "fail to reject" else "reject" end
    println(io, "Test summary:")
    label_len = length(conf_string) + 25 # 25 is length of 'outcome with  confidence:'
    println(io, "    outcome with $conf_string confidence: $outcome h_0")
    if tail == :both
        println(io, "    $(with_trailing_whitespace("two-sided p-value:", label_len)) $p")
    elseif tail == :left || tail == :right
        println(io, "    $(with_trailing_whitespace("one-sided p-value:", label_len)) $p")
    else
        println(io, "    $(with_trailing_whitespace("p-value:", label_len)) $p")
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
include("augmented_dickey_fuller.jl")
include("jarque_bera.jl")
end
