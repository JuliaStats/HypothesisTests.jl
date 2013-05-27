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
import Base.repl_show

export testname, pvalue, ci
abstract HypothesisTest

check_same_length(x::Vector, y::Vector) = if length(x) != length(y)
		error("Vectors must be the same length")
	end

# Basic function for finding a p-value given a distribution and tail
pvalue(dist::Distributions.Distribution, x::Number; tail=:both) = 
	if tail == :both
		2 * min(cdf(dist, x), ccdf(dist, x), 0.5)
	elseif tail == :left
		cdf(dist, x)
	elseif tail == :right
		ccdf(dist, x)
	else
		error("tail=$(tail) is invalid")
	end

# Repl pretty-print
function repl_show{T <: HypothesisTest}(io::IO, test::T)
	print(io, "$(test_name(T))\n\n")
	n = length(T.names)
	for i  = 1:n
		name = T.names[i]
		print(io, replace(string(name), "_", " "))
		print(io, " = $(getfield(test, name))")
		if i != n
			print(io, ", ")
		end
	end
end

include("circular.jl")
include("wilcoxon.jl")
include("t.jl")
end