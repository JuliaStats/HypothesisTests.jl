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

# TODO(cs): describe fisher in more detail (checks odds ratio)
# TODO(cs): confidence interval for odds ratio
# TODO(cs): Implement Blaker test
testname(::FisherExactTest) = "Fisher's exact test"
population_param_of_interest(x::FisherExactTest) = ("Odds ratio", 1, x.a/x.c/x.b*x.d) # parameter of interest: name, value under h0, point estimate

function show_params(io::IO, x::FisherExactTest, ident="")
	println(io, ident, "contingency table:")
	Base.print_matrix(io, [x.a x.b; x.c x.d], (typemax(Int), typemax(Int)), repeat(ident, 2))
	println(io)
end

# DOC: for tail=:both there exist multiple ``method``s for computing a pvalue and thus for the corresponding ci.
function pvalue(x::FisherExactTest; tail=:both, method=:central)
	if tail == :both && method != :central
		if method == :minlike
			pvalue_both_minlike(x)
		else
			error("method=$(method) is not implemented yet")
		end
	else
		pvalue(Hypergeometric(x.a + x.b, x.c + x.d, x.a + x.c), x.a, tail=tail)
	end
end

function pvalue_both_minlike(x::FisherExactTest)
	a, b, c, d = x.a, x.b, x.c, x.d

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
end

# confidence interval by inversion for central p-value
function ci_WIP(x::FisherExactTest, alpha::Float64=0.05; tail=:both, method=:central)
    check_alpha(alpha)

    if tail == :left
        (0.0, ci(x, alpha*2, method=:central)[2])
    elseif tail == :right
        (ci(x, alpha*2, method=:central)[1], inf)
    elseif tail == :both
	# needs noncentral hypergeometric distribution
        if method == :central
			error("method=$(method) is not implemented yet")
        else
            error("method=$(method) is not implemented yet")
        end
    else
        error("tail=$(tail) is invalid")
    end
end