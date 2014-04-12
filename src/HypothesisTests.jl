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
using Distributions

export testname, pvalue, ci
abstract HypothesisTest

check_same_length(x::Vector, y::Vector) = if length(x) != length(y)
		error("Vectors must be the same length")
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
		error("tail=$(tail) is invalid")
	end

pvalue(dist::DiscreteUnivariateDistribution, x::Number; tail=:both) = 
	if tail == :both
		min(2 * min(ccdf(dist, x-1), cdf(dist, x)), 1.0)
	elseif tail == :left
		cdf(dist, x)
	elseif tail == :right
		ccdf(dist, x-1)
	else
		error("tail=$(tail) is invalid")
	end

function check_alpha(alpha::Float64)
    if alpha <= 0 || alpha >= 0.5
        error("alpha $alpha not in range (0, 0.5)")
    end
end

# Pretty-print
function Base.show{T<:HypothesisTest}(io::IO, test::T)
	print(io, testname(test))
	fieldidx = find(Bool[t<:Number for t in T.types])
	if !isempty(fieldidx)
		print(io, "\n\n")
		lengths = [length(string(T.names[i])) for i in fieldidx]
		maxlen = maximum(lengths)

		for i = 1:length(fieldidx)
			name = T.names[fieldidx[i]]
			print(io, repeat(" ", maxlen-lengths[i]),
			          replace(string(name), "_", " "),
			          " = $(getfield(test, name))")
			if i != length(fieldidx)
				print(io, "\n")
			end
		end
	end

	has_pval = applicable(pvalue, test)
	has_ci = applicable(ci, test)
	if has_pval || has_ci
		println(io)
	end
	if has_pval
		print(io, "\nTwo-sided p-value:\n    p = $(pvalue(test))")
	end
	if has_ci
		confint = ci(test)
		print(io, "\n95% confidence interval:\n    $confint")
	end
end

include("binomial.jl")
include("circular.jl")
include("fisher.jl")
include("t.jl")
include("wilcoxon.jl")
end
