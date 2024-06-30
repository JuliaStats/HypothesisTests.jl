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
using Printf: @printf

import StatsAPI
using StatsAPI: HypothesisTest, confint, pvalue

export testname, pvalue, confint, dof, nobs

check_same_length(x::AbstractVector, y::AbstractVector) = if length(x) != length(y)
    throw(DimensionMismatch("Vectors must be the same length"))
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
    has_ci = applicable(confint, test)
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
        ci = map(x -> round.(x; sigdigits=4, base=10), confint(test))
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
include("kpss.jl")

end
