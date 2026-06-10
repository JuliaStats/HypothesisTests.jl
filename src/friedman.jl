export FriedmanTest, NemenyiTest, pvalues

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Compute within-row average ranks (tied values get the average rank).
# Returns an n × k matrix.
function _row_ranks(data::AbstractMatrix{<:Real})
	n, k = size(data)
	ranks = similar(data, Float64)
	for i in 1:n
		ranks[i, :] = tiedrank(view(data, i, :))
	end
	return ranks
end

# ---------------------------------------------------------------------------
# FriedmanTest
# ---------------------------------------------------------------------------

struct FriedmanTest <: HypothesisTest
	n::Int                      # number of blocks (datasets)
	k::Int                      # number of treatments (algorithms)
	ranks::Matrix{Float64}      # n × k matrix of within-block ranks
	avg_ranks::Vector{Float64}  # mean rank per treatment, length k
	chisq::Float64              # original Friedman χ²_F statistic
	FF::Float64                 # Iman-Davenport F statistic
	df1::Int                    # numerator df for F
	df2::Int                    # denominator df for F
end

"""
	FriedmanTest(data::AbstractMatrix{<:Real})

Perform the Friedman rank sum test on `data`, an `n × k` matrix where rows
are blocks (e.g. datasets) and columns are treatments (e.g. algorithms).

The null hypothesis is that all `k` treatments have equal distributions across
the `n` blocks (i.e. their expected average ranks are all equal).

Within each block the `k` values are replaced by their ranks 1 to k (ties
receive the average rank). The original Friedman statistic is

```math
\\chi^2_F = \\frac{12n}{k(k+1)} \\sum_{j=1}^{k}
			\\!\\left(\\bar{r}_j - \\frac{k+1}{2}\\right)^{\\!2}
```

distributed asymptotically as ``\\chi^2(k-1)``.  Demšar (2006) recommends
the Iman-Davenport F approximation, which has better Type-I error control:

```math
F_F = \\frac{(n-1)\\,\\chi^2_F}{n(k-1) - \\chi^2_F}
	  \\sim F\\bigl(k-1,\\,(n-1)(k-1)\\bigr)
```

`pvalue` returns the F-based p-value by default.
Pass `method = :chisq` for the chi-squared approximation.

Implements: [`pvalue`](@ref)

# References

  * [Demšar (2006). Statistical Comparisons of Classifiers over Multiple Data
	Sets](@cite demsar2006)
  * [Iman and Davenport (1980). Approximations of the critical region of the
	Friedman statistic](@cite iman1980)

# External links

  * [Friedman test on Wikipedia](https://en.wikipedia.org/wiki/Friedman_test)
"""
function FriedmanTest(data::AbstractMatrix{<:Real})
	n, k = size(data)
	n >= 2 || throw(ArgumentError("Need at least 2 blocks (rows), got $n"))
	k >= 2 || throw(ArgumentError("Need at least 2 treatments (columns), got $k"))

	ranks     = _row_ranks(data)
	avg_ranks = vec(mean(ranks; dims = 1))

	# Friedman χ²_F  (equivalent formulation via sum of squared avg ranks)
	chisq = 12n / (k * (k + 1)) * sum((avg_ranks .- (k + 1) / 2) .^ 2)

	# Iman-Davenport F  (Eq. 2 in Demšar 2006)
	denom = n * (k - 1) - chisq
	FF    = abs(denom) < eps(Float64) ? Inf : (n - 1) * chisq / denom

	return FriedmanTest(n, k, ranks, avg_ranks, chisq, FF, k - 1, (n - 1) * (k - 1))
end

testname(::FriedmanTest) = "Friedman rank sum test"
population_param_of_interest(::FriedmanTest) =
	("Average ranks of treatments", "all equal", NaN)
default_tail(::FriedmanTest) = :right

function StatsAPI.pvalue(x::FriedmanTest; method::Symbol = :f)
	if method === :f
		return pvalue(FDist(x.df1, x.df2), x.FF; tail = :right)
	elseif method === :chisq
		return pvalue(Chisq(x.df1), x.chisq; tail = :right)
	else
		throw(ArgumentError("method must be :f or :chisq, got :$method"))
	end
end

function show_params(io::IO, x::FriedmanTest, ident)
	println(io, ident, "number of blocks (n):       ", x.n)
	println(io, ident, "number of treatments (k):   ", x.k)
	print(io, ident, "average ranks:              ")
	show(io, x.avg_ranks);
	println(io)
	println(io, ident, "χ²_F statistic:             ", x.chisq)
	println(io, ident, "Iman-Davenport F statistic: ", x.FF)
	println(io, ident, "F degrees of freedom:       (", x.df1, ", ", x.df2, ")")
end


# ---------------------------------------------------------------------------
# NemenyiTest
# ---------------------------------------------------------------------------

# q_α critical values from Studentized range distribution / √2
# Tabulated in Demšar (2006) Table 5 for k = 2…10, α ∈ {0.10, 0.05, 0.01}.
# Source: Zar (1999) "Biostatistical Analysis", Table B.5.
const _NEMENYI_Q = Dict{Int, NTuple{3, Float64}}(
	# q_{α,k} = StudentizedRange(k, ∞) / √2, matching Demšar (2006) Table 5.
	# Columns: α = 0.10, 0.05, 0.01
	2  => (1.6449, 1.9600, 2.5758),
	3  => (2.0523, 2.3437, 2.9135),
	4  => (2.2913, 2.5690, 3.1133),
	5  => (2.4595, 2.7278, 3.2547),
	6  => (2.5885, 2.8497, 3.3637),
	7  => (2.6927, 2.9483, 3.4522),
	8  => (2.7799, 3.0309, 3.5265),
	9  => (2.8546, 3.1017, 3.5903),
	10 => (2.9199, 3.1637, 3.6463),
)

struct NemenyiTest <: HypothesisTest
	n::Int                      # number of blocks
	k::Int                      # number of treatments
	avg_ranks::Vector{Float64}  # mean rank per treatment
	cd::Float64                 # critical difference at chosen α
	pvalues::Matrix{Float64}    # k × k symmetric p-value matrix
	alpha::Float64              # significance level
end

"""
	NemenyiTest(ft::FriedmanTest; alpha::Real = 0.05)

Perform the Nemenyi all-pairs post-hoc test following a significant
[`FriedmanTest`](@ref), as recommended by Demšar (2006).

Two treatments `i` and `j` differ significantly when the absolute difference
in their average ranks exceeds the critical difference

```math
CD = q_\\alpha \\sqrt{\\frac{k(k+1)}{6n}}
```

where ``q_\\alpha`` is tabulated from the Studentized range distribution divided
by ``\\sqrt{2}`` (Demšar 2006, Table 5; computed as
``\\text{StudentizedRange}(k, \\infty).\\text{ppf}(1-\\alpha) / \\sqrt{2}``; tabulated
for ``k = 2, \\ldots, 10`` and ``\\alpha \\in \\{0.10, 0.05, 0.01\\}``).

Pairwise p-values are approximated via the Bonferroni-adjusted normal
approximation to the Studentized range distribution at infinite degrees of
freedom, which is the same approximation used by Demšar (2006):

```math
p_{ij} = \\min\\!\\left(1,\\; k(k-1) \\cdot 2\\Phi\\!\\left(
		  -\\frac{|\\bar{r}_i - \\bar{r}_j|}{\\sqrt{k(k+1)/6n}}
		  \\right)\\right)
```

Implements: [`pvalue`](@ref)

!!! note
	Tabulated ``q_\\alpha`` values are available for ``k \\leq 10``.  For larger
	``k`` the Bonferroni-adjusted normal quantile is used as a conservative
	approximation.

# References

  * [Demšar (2006). Statistical Comparisons of Classifiers over Multiple Data
	Sets](@cite demsar2006)
  * [Nemenyi (1963). Distribution-free multiple comparisons. Ph.D. thesis,
	Princeton University.](@cite nemenyi1963)

# External links

  * [Nemenyi test on Wikipedia](https://en.wikipedia.org/wiki/Nemenyi_test)
"""
function NemenyiTest(ft::FriedmanTest; alpha::Real = 0.05)
	0 < alpha < 1 || throw(ArgumentError(lazy"alpha must be in (0, 1), got $alpha"))

	n, k      = ft.n, ft.k
	avg_ranks = ft.avg_ranks
	m         = k * (k - 1)          # total number of pairwise comparisons × 2
	se        = sqrt(m / (6 * n))

	# Bonferroni-adjusted two-sided p-values using the Normal approximation
	# to the Studentized range distribution at ∞ degrees of freedom.
	# This matches the approximation implicit in Demšar (2006).
	pvalues = zeros(Float64, k, k)
	for i in 1:k
		for j in (i+1):k
			z = abs(avg_ranks[i] - avg_ranks[j]) / se
			p = min(1.0, m * 2 * cdf(Normal(), -z))
			pvalues[i, j] = p
			pvalues[j, i] = p
		end
	end

	# Critical difference from tabulated q_α (or Bonferroni normal fallback)
	q = if haskey(_NEMENYI_Q, k)
		vals = _NEMENYI_Q[k]
		alpha <= 0.01 + eps() ? vals[3] :
		alpha <= 0.05 + eps() ? vals[2] :
		vals[1]
	else
		# Conservative fallback: Bonferroni normal quantile
		quantile(Normal(), 1 - alpha / m)
	end

	return NemenyiTest(n, k, avg_ranks, q * se, pvalues, Float64(alpha))
end

testname(::NemenyiTest) = "Nemenyi all-pairs post-hoc test"
population_param_of_interest(::NemenyiTest) =
	("Pairwise average rank differences", "all zero", NaN)
default_tail(::NemenyiTest) = :both

"""
	pvalue(x::NemenyiTest)

Return the minimum pairwise Bonferroni-adjusted p-value across all treatment pairs.
"""
StatsAPI.pvalue(x::NemenyiTest) =
	minimum(x.pvalues[i, j] for i in 1:x.k for j in 1:x.k if i != j)

"""
	pvalue(x::NemenyiTest, i::Int, j::Int)

Return the p-value for the pairwise comparison of treatments `i` and `j`.
"""
StatsAPI.pvalue(x::NemenyiTest, i::Int, j::Int) = x.pvalues[i, j]

"""
	pvalues(x::NemenyiTest)

Return the full ``k × k`` matrix of pairwise Bonferroni-adjusted p-values.
Diagonal entries are 0; the matrix is symmetric.
"""
pvalues(x::NemenyiTest) = x.pvalues

function show_params(io::IO, x::NemenyiTest, ident)
	println(io, ident, "number of blocks (n):      ", x.n)
	println(io, ident, "number of treatments (k):  ", x.k)
	print(io, ident, "average ranks:             ")
	show(io, x.avg_ranks);
	println(io)
	println(io, ident, "significance level (α):    ", x.alpha)
	println(io, ident, "critical difference (CD):  ", x.cd)
	println(io, ident, "pairwise p-value matrix (Bonferroni-adjusted):")
	for i in 1:x.k
		print(io, ident, "  ")
		for j in 1:x.k
			@printf(io, "  %7.4f", x.pvalues[i, j])
		end
		println(io)
	end
end
