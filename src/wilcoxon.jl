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

export MannWhitneyUTest, ExactMannWhitneyUTest, ApproximateMannWhitneyUTest,
	SignedRankTest, ExactSignedRankTest, ApproximateSignedRankTest

## EXACT MANN-WHITNEY U TEST

abstract MannWhitneyUTest <: HypothesisTest
type ExactMannWhitneyUTest <: MannWhitneyUTest
	U::Float64
	p_value::Float64
end
test_name(::Type{ExactMannWhitneyUTest}) = "Exact Mann-Whitney U test"

# Tied rank from Base, modified to compute the adjustment for ties
function tiedrank_adj(v::AbstractArray)
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

    (ord, tieadj)
end

# Enumerate all possible Mann-Whitney U results for a given vector, determining left-
# and right-tailed p values
function mwu_enumerate{S <: Real}(nx::Int, ny::Int, U::Float64, ranks::Vector{S})
	# Get the other U if inverted by mwu_stats
	n = min(nx, ny)
	if ny > nx
		U = nx*ny - U
	end
	le = 0
	gr = 0
	tot = 0
	k = n*(n+1)/2
	for comb in @task combinations(ranks, n)
		Up = sum(comb) - k
		tot += 1
		le += Up <= U
		gr += Up >= U
	end
	(le/tot, gr/tot)
end

p_value{S <: Real}(::Type{ExactMannWhitneyUTest}, U::Real, ranks::Vector{S}, tieadj::Int, nx::Int, ny::Int) =
	if tieadj == 0
		# Compute exact p-value using method from Rmath, which is fast but cannot account for ties
		if U < nx * ny / 2
			2 * pwilcox(U, nx, ny, true)
		else
			2 * pwilcox(U - 1, nx, ny, false)
		end
	else
		# Compute exact p-value by enumerating all possible ranks in the tied data
		min(1, 2 * min(mwu_enumerate(nx, ny, U, ranks)))
	end
left_p_value{S <: Real}(::Type{ExactMannWhitneyUTest}, U::Real, ranks::Vector{S}, tieadj::Int, nx::Int, ny::Int) =
	tieadj == 0 ? pwilcox(U, nx, ny, true) : mwu_enumerate(nx, ny, U, ranks)[1]
right_p_value{S <: Real}(::Type{ExactMannWhitneyUTest}, U::Real, ranks::Vector{S}, tieadj::Int, nx::Int, ny::Int) =
	tieadj == 0 ? pwilcox(U - 1, nx, ny, false) : mwu_enumerate(nx, ny, U, ranks)[2]

## APPROXIMATE MANN-WHITNEY U TEST

type ApproximateMannWhitneyUTest <: MannWhitneyUTest
	U::Float64
	p_value::Float64
end
test_name(::Type{ApproximateMannWhitneyUTest}) = "Approximate Mann-Whitney U test"

# Get mean and sigma for null distribution of approximate Mann-Whitney U test (without continuity correction)
mwu_z(U::Real, tieadj::Int, nx::Int, ny::Int) =
	(U - nx * ny / 2, sqrt((nx * ny * (nx + ny + 1 - tieadj / ((nx + ny) * (nx + ny - 1)))) / 12))

let d = Normal()
	function p_value{S <: Real}(::Type{ApproximateMannWhitneyUTest}, U::Real, ranks::Vector{S}, tieadj::Int, nx::Int, ny::Int)
		(mu, sigma) = mwu_z(U, tieadj, nx, ny)
		2 * ccdf(d, abs(mu - 0.5 * sign(mu))/sigma)
	end
	function left_p_value{S <: Real}(::Type{ApproximateMannWhitneyUTest}, U::Real, ranks::Vector{S}, tieadj::Int, nx::Int, ny::Int)
		(mu, sigma) = mwu_z(U, tieadj, nx, ny)
		cdf(d, (mu + 0.5)/sigma)
	end
	function right_p_value{S <: Real}(::Type{ApproximateMannWhitneyUTest}, U::Real, ranks::Vector{S}, tieadj::Int, nx::Int, ny::Int)
		(mu, sigma) = mwu_z(U, tieadj, nx, ny)
		ccdf(d, (mu - 0.5)/sigma)
	end
end

## COMMON MANN-WHITNEY U

# Get U, ranks, and tie adjustment for Mann-Whitney U test
function mwu_stats{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T})
	nx = length(x)
	ny = length(y)
	if nx <= ny
		(ranks, tieadj) = tiedrank_adj([x, y])
		U = sum(ranks[1:nx]) - nx*(nx+1)/2
	else
		(ranks, tieadj) = tiedrank_adj([y, x])
		U = nx*ny - sum(ranks[1:ny]) + ny*(ny+1)/2
	end
	(U, ranks, tieadj, nx, ny)
end

# Constructors
for t in (:ExactMannWhitneyUTest, :ApproximateMannWhitneyUTest)
	@eval begin
		function $(t){S <: Real, T <: Real}(x::Vector{S}, y::Vector{T})
			(U, ranks, tieadj, nx, ny) = mwu_stats(x, y)
			p = p_value($(t), U, ranks, tieadj, nx, ny)
			$(t)(U, p)
		end
	end
end

# Test statistic and p-values
test_statistic{S <: Real, T <: Real, U <: MannWhitneyUTest}(::Type{U}, x::Vector{S}, y::Vector{T}) =
	mwu_stats(x, y)[1]
for fn in (:p_value, :left_p_value, :right_p_value)
	@eval begin
		$(fn){S <: Real, T <: Real, U <: MannWhitneyUTest}(::Type{U}, x::Vector{S}, y::Vector{T}) =
			$(fn)(U, mwu_stats(x, y)...)
	end
end

# Automatic exact/normal selection
for fn in (:p_value, :left_p_value, :right_p_value, :test_statistic)
	@eval begin
		function $(fn){S <: Real, T <: Real}(::Type{MannWhitneyUTest}, x::Vector{S}, y::Vector{T})
			(U, ranks, tieadj, nx, ny) = mwu_stats(x, y)
			if nx + ny <= 10 || (nx + ny <= 50 && tieadj == 0)
				$(fn)(ExactMannWhitneyUTest, U, ranks, tieadj, nx, ny)
			else
				$(fn)(ApproximateMannWhitneyUTest, U, ranks, tieadj, nx, ny)
			end
		end
	end
end

## EXACT WILCOXON SIGNED RANK TEST

abstract SignedRankTest <: HypothesisTest
type ExactSignedRankTest <: SignedRankTest
	W::Float64
	p_value::Float64
end
test_name(::Type{ExactSignedRankTest}) = "Exact Wilcoxon signed rank test"

# Enumerate all possible Mann-Whitney U results for a given vector, determining left-
# and right-tailed p values
function signrank_enumerate{T <: Real}(W::Real, ranks::Vector{T})
	le = 0
	gr = 0
	n = length(ranks)
	tot = 2^n
	for i = 0:tot-1
		# Interpret bits of i as signs to generate wp for all possible sign combinations
		Wp = 0
		x = i
		j = 1
		while x != 0
			Wp += (x & 1)*ranks[j]
			j += 1
			x >>= 1
		end
		le += Wp <= W
		gr += Wp >= W
	end
	(le/tot, gr/tot)
end

p_value{T <: Real}(::Type{ExactSignedRankTest}, W::Real, ranks::Vector{T}, tieadj::Int) = 
	if length(ranks) == 0
		1
	elseif tieadj == 0
		# Compute exact p-value using method from Rmath, which is fast but cannot account for ties
		n = length(ranks)
		if W <= n * (n + 1)/4
			p = 2 * psignrank(W, n, true)
		else
			p = 2 * psignrank(W - 1, n, false)
		end
	else
		# Compute exact p-value by enumerating all possible ranks in the tied data
		min(2 * min(signrank_enumerate(W, ranks)...), 1)
	end
left_p_value{S <: Real}(::Type{ExactSignedRankTest}, W::Real, ranks::Vector{S}, tieadj::Int) =
	length(ranks) == 0 ? 1 : tieadj == 0 ? psignrank(W, length(ranks), true) : signrank_enumerate(W, ranks)[1]
right_p_value{S <: Real}(::Type{ExactSignedRankTest}, W::Real, ranks::Vector{S}, tieadj::Int) =
	length(ranks) == 0 ? 1 : tieadj == 0 ? psignrank(W - 1, length(ranks), false) : signrank_enumerate(W, ranks)[2]

## APPROXIMATE SIGNED RANK TEST

type ApproximateSignedRankTest <: SignedRankTest
	W::Float64
	p_value::Float64
end
test_name(::Type{ApproximateSignedRankTest}) = "Approximate Wilcoxon signed rank test"

# Get mean and sigma for null distribution of approximate signed rank test (without continuity correction)
signrank_z(W::Real, n::Int, tieadj::Int) =
	(W - n * (n + 1)/4, sqrt(n * (n + 1) * (2 * n + 1) / 24 - tieadj / 48))

let d = Normal()
	function p_value{S <: Real}(::Type{ApproximateSignedRankTest}, W::Real, ranks::Vector{S}, tieadj::Int)
		if length(ranks) == 0
			return 1
		end
		(mu, sigma) = signrank_z(W, length(ranks), tieadj)
		2 * ccdf(d, abs(mu - 0.5 * sign(mu))/sigma)
	end
	function left_p_value{S <: Real}(::Type{ApproximateSignedRankTest}, W::Real, ranks::Vector{S}, tieadj::Int)
		if length(ranks) == 0
			return 1
		end
		(mu, sigma) = signrank_z(W, length(ranks), tieadj)
		cdf(d, (mu + 0.5)/sigma)
	end
	function right_p_value{S <: Real}(::Type{ApproximateSignedRankTest}, W::Real, ranks::Vector{S}, tieadj::Int)
		if length(ranks) == 0
			return 1
		end
		(mu, sigma) = signrank_z(W, length(ranks), tieadj)
		ccdf(d, (mu - 0.5)/sigma)
	end
end

## COMMON SIGNED RANK

# Get W and absolute ranks for signed rank test
function signrank_stats{S <: Real}(x::Vector{S})
	nonzero_x = x[x .!= 0]
	(ranks, tieadj) = tiedrank_adj(abs(nonzero_x))
	W = 0.0
	for i = 1:length(nonzero_x)
		if nonzero_x[i] > 0
			W += ranks[i]
		end
	end
	(W, ranks, tieadj)
end

# Constructors
for t in (:ExactSignedRankTest, :ApproximateSignedRankTest)
	@eval begin
		function $(t){S <: Real}(x::Vector{S})
			(W, ranks, tieadj) = signrank_stats(x)
			p = p_value($(t), W, ranks, tieadj)
			$(t)(W, p)
		end

		$(t){S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}) = $(t)(x - y)
	end
end

# Test statistic and p-values
test_statistic{S <: Real, U <: SignedRankTest}(::Type{U}, x::Vector{S}) = signrank_stats(x)[1]
test_statistic{S <: Real, T <: Real, U <: SignedRankTest}(::Type{U}, x::Vector{S}, y::Vector{T}) =
	signrank_stats(x - y)[1]
for fn in (:p_value, :left_p_value, :right_p_value)
	@eval begin
		$(fn){S <: Real, U <: SignedRankTest}(::Type{U}, x::Vector{S}) =
			$(fn)(U, signrank_stats(x)...)
		$(fn){S <: Real, T <: Real, U <: SignedRankTest}(::Type{U}, x::Vector{S}, y::Vector{T}) =
			$(fn)(U, signrank_stats(x - y)...)
	end
end

# Automatic exact/normal selection
for fn in (:p_value, :left_p_value, :right_p_value, :test_statistic)
	@eval begin
		function $(fn){S <: Real}(::Type{SignedRankTest}, W::Real, ranks::Vector{S}, tieadj::Int)
			n = length(ranks)
			if n <= 15 || (n <= 50 && tieadj == 0)
				$(fn)(ExactSignedRankTest, W, ranks, tieadj)
			else
				$(fn)(ApproximateSignedRankTest, W, ranks, tieadj)
			end
		end
	end
end