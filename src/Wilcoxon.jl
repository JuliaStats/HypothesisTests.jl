# Wilcoxon.jl
# Wilcoxon rank sum (Mann-Whitney U) and signed rank tests in Julia
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

load("Distributions")

module Wilcoxon
include("/usr/local/julia/extras/Rmath.jl")
using Distributions

export mannwhitneyu, TWO_TAILED, LEFT_TAILED, RIGHT_TAILED

const TWO_TAILED = 0
const LEFT_TAILED = 1
const RIGHT_TAILED = 2

# order (aka, rank), resolving ties using the mean rank and computing adjustment
function tiedrank(v::AbstractArray)
    n     = length(v)
    place = order(v)
    ord   = Array(Float64, n)
    tieadj = 0

    i = 1
    while i <= n
        j = i
        while j + 1 <= n && v[place[i]] == v[place[j + 1]]
            j += 1
        end

        if j > i
            t = j - i + 1
            m = sum(i:j) / t
            tieadj += t^3 - t
            for k = i:j
                ord[place[k]] = m
            end
        else
            ord[place[i]] = i
        end

        i = j + 1
    end

    return (ord, tieadj)
end

# Mann-Whitney U test
function mannwhitneyu{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}, tail::Int)
	nx = length(x)
	ny = length(y)
	n_total = nx + ny
	if nx < ny
		(ranks, tieadj) = tiedrank([x, y])
		n = nx
	else
		(ranks, tieadj) = tiedrank([y, x])
		n = ny
		if tail == LEFT_TAILED
			tail = RIGHT_TAILED
		elseif tail == RIGHT_TAILED
			tail = LEFT_TAILED
		end
	end

	U = sum(ranks[1:n]) - n*(n+1)/2
	
	if n_total <= 50 && tieadj == 0
		# Compute exact p-value using method from Rmath, which is fast but cannot account
		# for ties in the data
		if tail == LEFT_TAILED
			p = pwilcox(U, nx, ny, true)
		elseif tail == RIGHT_TAILED
			p = pwilcox(U - 1, nx, ny, false)
		else
			if U < nx * ny / 2
				p = 2 * pwilcox(U, nx, ny, true)
			else
				p = 2 * pwilcox(U - 1, nx, ny, false)
			end
		end
	elseif n_total <= 10
		# Compute exact p-value by enumerating all possible ranks in the tied data
		le = 0
		gr = 0
		tot = 0
		k = n*(n+1)/2
		for comb in @task combinations(ranks, n)
			Up = sum(comb[1:n]) - k
			tot += 1
			le += Up <= U
			gr += Up >= U
		end
		if tail == LEFT_TAILED
			p = le/tot
		elseif tail == RIGHT_TAILED
			p = gr/tot
		else
			p = 2 * min([le, gr]/tot)
		end
	else
		# Compute approximate p-value
		d = Normal(nx * ny / 2,
			sqrt((nx * ny * (n_total + 1 - tieadj / (n_total * (n_total - 1)))) / 12))
		if tail == LEFT_TAILED
			p = cdf(d, U + 1/2)
		elseif tail == RIGHT_TAILED
			p = ccdf(d, U - 1/2)
		else
			p = 2 * (U < nx * ny / 2 ? cdf(d, U + 1/2) : ccdf(d, U - 1/2))
		end
	end
	return (p, U)
end

# Wilcoxon signed rank test
function signrank{S <: Real}(x::Vector{S}, tail::Int)
	n = length(x)

	(ranks, tieadj) = tiedrank(x)
	absranks = abs(ranks)
	w = sum(absranks .* sign(ranks))
	
	if n <= 50 && tieadj == 0
		# Compute exact p-value using method from Rmath, which is fast but cannot account
		# for ties in the data
		if tail == LEFT_TAILED
			p = psignrank(w, n, true)
		elseif tail == RIGHT_TAILED
			p = psignrank(w - 1, n, false)
		else
			if w < 0
				p = 2 * psignrank(w, n, true)
			else
				p = 2 * psignrank(w - 1, n, false)
			end
		end
	elseif n <= 15
		# Compute exact p-value by enumerating all possible signs in the tied data
		le = 0
		gr = 0
		tot = 2^n
		for i = 0:tot-1
			# Interpret bits of i as signs to generate wp for all possible sign
			# combinations
			x = i
			wp = 0
			for j = 1:n
				wp += ((x & 1 == 0) * 2 - 1) * absranks[j]
				x >>= 1
			end
			le += wp <= w
			gr += wp >= w
		end
		if tail == LEFT_TAILED
			p = le/tot
		elseif tail == RIGHT_TAILED
			p = gr/tot
		else
			p = 2 * min([le, gr]/tot)
		end
	else
		# Compute approximate p-value
		d = Normal(n * (n + 1)/4,
			sqrt(n * (n + 1) * (2 * n + 1) / 24 - tieadj / 48))
		if tail == LEFT_TAILED
			p = cdf(d, w + 1/2)
		elseif tail == RIGHT_TAILED
			p = cdf(d, -w - 1/2)
		else
			p = 2 * (w < n * (n + 1)/4 ? cdf(d, w + 1/2) : ccdf(d, w - 1/2))
		end
	end
	return (p, w)
end
signrank{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}) =
	signrank(x - y, TWO_TAILED)

end