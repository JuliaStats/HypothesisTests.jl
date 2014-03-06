# fisher.jl
# Fisher's exact test
#
# Copyright (C) 2013   Simon Kornblith
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

export FisherExactTest

immutable FisherExactTest <: HypothesisTest
	# Format:
	#    X1  X2
	# Y1  a  b
	# Y2  c  d
	a::Int
	b::Int
	c::Int
	d::Int
end

testname(::FisherExactTest) = "Fisher's exact test"
function pvalue(x::FisherExactTest; tail=:both)
	a, b, c, d = x.a, x.b, x.c, x.d

	if tail == :both
		if a + c > b + d
			a, b, c, d = b, a, d, c
		end
		if a/c > b/d
			a, b, c, d = c, d, a, b
		end
		dist = Hypergeometric(a+b, c+d, a+c)

		p = pdf(dist, a)
		v = nextfloat(p)
		if a != 0
			p += cdf(dist, a-1)
		end

		# Add p-values of all tables in other tail equally or less probable
		for i = a+c:-1:a+1
			curp = pdf(dist, i)
			if curp > v
				break
			end
			p += curp
		end
		p
	else
		pvalue(Hypergeometric(a+b, c+d, a+c), a, tail=tail)
	end
end

function Base.show(io::IO, x::FisherExactTest)
	a, b, c, d = x.a, x.b, x.c, x.d
	print(io, testname(x), "\n\nContingency table:\n")
	Base.print_matrix(io, [a b; c d], typemax(Int), typemax(Int), "    ")
	print(io, """\n\n
Two-sided p-value:
    p = $(pvalue(x))

Odds ratio:
    $(a/c/b*d)""")
end
