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
    
    # utilities for pretty-printing
    conf_string = string(floor((1 - alpha(test)) * 100, 6)) # limit to 6 decimals in %
    prettify_detail(label::String, value::Any, len::Int) = # len is max length of label
        "    " * label * " "^max(len - length(label), 0) * string(value)

    # population details
    has_ci = applicable(StatsBase.confint, test)
    (param_name, param_under_h0, param_estimate) = population_param_of_interest(test)
    println(io, "Population details:")
    println(io, prettify_detail("parameter of interest:", param_name, 32))
    println(io, prettify_detail("value under h_0:", param_under_h0, 32))
    println(io, prettify_detail("point estimate:", param_estimate, 32))
    if has_ci
        println(io, prettify_detail(conf_string*"% confidence interval:", StatsBase.confint(test), 32))
    end
    println(io)

    # test summary
    tail = HypothesisTests.tail(test)
    p = pvalue(test) # obeys value of HypothesisTests.tail(test) if applicable
    outcome = if p > alpha(test) "fail to reject" else "reject" end
    tailvalue =
        if     tail == :both "two-sided p-value:"
        elseif tail == :left || tail == :right "one-sided p-value ($(string(tail)) tail):"
        else   "p-value:" end
    println(io, "Test summary:")
    println(io, prettify_detail("outcome with "*conf_string*"% confidence:", outcome*" h_0", 36))
    println(io, prettify_detail(tailvalue, p, 36))
    println(io)

    # further details
    println(io, "Details:")
    show_params(io, test, "    ")
end

# parameter of interest: name, value under h0, point estimate
population_param_of_interest{T<:HypothesisTest}(test::T) = ("not implemented yet", NaN, NaN)

# is the test one- or two-sided?
tail(test::HypothesisTest)  = :undefined    # overloaded for defaults or field access
alpha(test::HypothesisTest) = 0.05

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
