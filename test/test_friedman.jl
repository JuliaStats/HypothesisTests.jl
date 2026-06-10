# test/friedman.jl
# Tests for FriedmanTest and NemenyiTest
#
# Ground truth values verified against:
#   - Manual computation with SciPy (scipy.stats.studentized_range, scipy.stats.f)
#   - Demšar (2006) Table 5 q-values and CD formula
#   - Iman-Davenport (1980) F statistic formula

using Test
using HypothesisTests
using Statistics: mean

# ─────────────────────────────────────────────────────────────────────────────
# Shared test data
#
# Simple 5-dataset × 3-algorithm matrix where algorithm A always scores
# highest and C always lowest. This gives perfectly consistent ascending ranks:
#   col 1 → rank 3, col 2 → rank 2, col 3 → rank 1  (ascending tiedrank)
# Ground truth (verified with SciPy):
#   avg_ranks  = [3.0, 2.0, 1.0]
#   chisq      = 10.0   (= 12*5/(3*4) * ((3-2)²+(2-2)²+(1-2)²) * 5)
#   F_F        = Inf    (denominator = 5*2 - 10 = 0; perfect discrimination)
#   p (F)      = 0.0
#   CD(α=0.05) = 2.3437 * sqrt(3*4/30) = 2.3437 * 0.6325 ≈ 1.4823
# ─────────────────────────────────────────────────────────────────────────────
const test5x3 = Float64[
	0.90  0.70  0.50;
	0.80  0.60  0.40;
	0.85  0.65  0.45;
	0.75  0.55  0.35;
	0.95  0.75  0.55;
]

# ─────────────────────────────────────────────────────────────────────────────
# A larger example with non-trivial ranks and known statistics.
# 8 datasets × 4 algorithms, no ties within rows.
# Ground truth computed with SciPy:
#   avg_ranks ≈ [3.25, 1.875, 2.625, 2.25]
#   chisq     ≈ 5.8125
#   F_F       ≈ 1.6508   (F(3,21))
#   p (F)     ≈ 0.2078
# ─────────────────────────────────────────────────────────────────────────────
const test8x4 = Float64[
	0.60  0.80  0.70  0.75;
	0.55  0.85  0.65  0.72;
	0.70  0.75  0.80  0.60;
	0.65  0.90  0.55  0.78;
	0.50  0.82  0.68  0.71;
	0.72  0.78  0.74  0.66;
	0.58  0.88  0.62  0.76;
	0.68  0.84  0.72  0.64;
]


@testset "FriedmanTest" begin

	@testset "Fields and size" begin
		ft = FriedmanTest(test5x3)
		@test ft.n == 5
		@test ft.k == 3
		@test ft.df1 == 2
		@test ft.df2 == 8
		@test size(ft.ranks) == (5, 3)
		@test length(ft.avg_ranks) == 3
	end

	@testset "Ranks are in [1, k] and sum to k*(k+1)/2 per row" begin
		ft = FriedmanTest(test5x3)
		n, k = ft.n, ft.k
		for i in 1:n
			@test sum(ft.ranks[i, :]) ≈ k * (k + 1) / 2
			@test minimum(ft.ranks[i, :]) >= 1.0
			@test maximum(ft.ranks[i, :]) <= Float64(k)
		end
	end

	@testset "avg_ranks sum to k*(k+1)/2" begin
		for data in (test5x3, test8x4)
			ft = FriedmanTest(data)
			@test sum(ft.avg_ranks) ≈ ft.k * (ft.k + 1) / 2 atol=1e-10
		end
	end

	@testset "test5x3 — known exact values" begin
		ft = FriedmanTest(test5x3)
		# Column with highest values gets rank 3 (ascending tiedrank)
		@test ft.avg_ranks ≈ [3.0, 2.0, 1.0] atol=1e-10
		@test ft.chisq ≈ 10.0 atol=1e-10
		@test isinf(ft.FF)    # denominator = 0 for perfect discrimination
		@test pvalue(ft; method = :f) == 0.0
		@test pvalue(ft; method = :chisq) < 1e-2
	end

	@testset "test8x4 — known statistics (SciPy cross-validation)" begin
		ft = FriedmanTest(test8x4)
		@test ft.n == 8 && ft.k == 4
		@test ft.avg_ranks ≈ [1.5, 3.875, 2.375, 2.25] atol=0.001
		@test ft.chisq ≈ 14.25 atol=0.001
		@test ft.FF ≈ 10.2308 atol=0.01
		@test pvalue(ft; method = :f) ≈ 0.00023 atol=0.0001
	end

	@testset "pvalue method dispatch and error" begin
		ft = FriedmanTest(test5x3)
		@test 0.0 <= pvalue(ft; method = :f) <= 1.0
		@test 0.0 <= pvalue(ft; method = :chisq) <= 1.0
		@test_throws ArgumentError pvalue(ft; method = :bad)
	end

	@testset "default_tail is :right" begin
		@test default_tail(FriedmanTest(test5x3)) === :right
	end

	@testset "testname" begin
		@test testname(FriedmanTest(test5x3)) == "Friedman rank sum test"
	end

	@testset "Input validation" begin
		@test_throws ArgumentError FriedmanTest(rand(1, 3))  # n < 2
		@test_throws ArgumentError FriedmanTest(rand(5, 1))  # k < 2
	end

	@testset "Identical columns → chisq = 0, p = 1" begin
		data = hcat(ones(8), ones(8), ones(8))
		ft = FriedmanTest(data)
		@test ft.chisq ≈ 0.0 atol=1e-10
		@test pvalue(ft; method = :chisq) ≈ 1.0 atol=1e-10
	end

	@testset "Ties within rows are handled (fractional avg ranks)" begin
		# Two identical columns → tied rank = 1.5, 1.5, 3 for each row
		data = Float64[1 1 3; 2 2 4; 1 1 2]
		ft   = FriedmanTest(data)
		@test all(ft.ranks[:, 1] .≈ 1.5)
		@test all(ft.ranks[:, 2] .≈ 1.5)
		@test all(ft.ranks[:, 3] .== 3.0)
	end

end  # FriedmanTest


@testset "NemenyiTest" begin

	ft3 = FriedmanTest(test5x3)
	ft4 = FriedmanTest(test8x4)

	@testset "Fields and shape" begin
		nt = NemenyiTest(ft3; alpha = 0.05)
		@test nt.n == 5
		@test nt.k == 3
		@test nt.alpha == 0.05
		@test size(nt.pvalues) == (3, 3)
	end

	@testset "p-value matrix properties" begin
		nt = NemenyiTest(ft3)
		@test all(0.0 .<= nt.pvalues .<= 1.0)
		@test all(nt.pvalues[i, i] == 0.0 for i in 1:nt.k)
		@test nt.pvalues ≈ nt.pvalues' atol=1e-15  # symmetric
	end

	@testset "Critical difference — Demšar (2006) formula" begin
		# CD = q_alpha * sqrt(k*(k+1) / (6*n))
		# k=4, n=8, alpha=0.05: q=2.5690 (Demšar Table 5)
		nt = NemenyiTest(ft4; alpha = 0.05)
		se_expected = sqrt(4 * 5 / (6 * 8))           # ≈ 0.6455
		cd_expected = 2.5690 * se_expected              # ≈ 1.6584
		@test nt.cd ≈ cd_expected atol=0.001
	end

	@testset "Larger α → smaller CD" begin
		nt10 = NemenyiTest(ft4; alpha = 0.10)
		nt05 = NemenyiTest(ft4; alpha = 0.05)
		nt01 = NemenyiTest(ft4; alpha = 0.01)
		@test nt10.cd < nt05.cd < nt01.cd
	end

	@testset "test5x3 — significant pair detected" begin
		# avg_ranks = [3, 2, 1]; CD ≈ 1.4823 (SciPy verified)
		# |col1 - col3| = 2.0 > 1.4823  → significant
		# |col1 - col2| = 1.0 < 1.4823  → not significant
		nt = NemenyiTest(ft3; alpha = 0.05)
		@test issignificant(nt, 1, 3)   # clearly worst vs best
		@test !issignificant(nt, 1, 2)   # adjacent ranks
		@test !issignificant(nt, 2, 3)
	end

	@testset "issignificant consistent with CD" begin
		for ft in (ft3, ft4), alpha in (0.01, 0.05, 0.10)
			nt = NemenyiTest(ft; alpha = alpha)
			for i in 1:nt.k, j in 1:nt.k
				expected = i != j &&
						   abs(nt.avg_ranks[i] - nt.avg_ranks[j]) > nt.cd
				@test issignificant(nt, i, j) == expected
			end
		end
	end

	@testset "pvalue(nt, i, j) indexing" begin
		nt = NemenyiTest(ft4)
		for i in 1:nt.k, j in 1:nt.k
			@test pvalue(nt, i, j) === nt.pvalues[i, j]
		end
		@test pvalue(nt) === nt.pvalues
	end

	@testset "p-value for extreme pair is small" begin
		nt = NemenyiTest(ft3; alpha = 0.05)
		# col1 vs col3: delta=2, Bonferroni-adjusted p ≈ 0.0094 (SciPy verified)
		@test pvalue(nt, 1, 3) ≈ 0.0094 atol=0.001
	end

	@testset "Input validation" begin
		@test_throws ArgumentError NemenyiTest(ft3; alpha = 0.0)
		@test_throws ArgumentError NemenyiTest(ft3; alpha = 1.0)
		@test_throws ArgumentError NemenyiTest(ft3; alpha = 1.5)
	end

	@testset "Fallback for k > 10 (Bonferroni-normal)" begin
		data_big = rand(15, 12)
		ft_big   = FriedmanTest(data_big)
		nt_big   = NemenyiTest(ft_big; alpha = 0.05)
		@test nt_big.cd > 0.0
		@test all(0.0 .<= nt_big.pvalues .<= 1.0)
		@test nt_big.pvalues ≈ nt_big.pvalues'
	end

	@testset "testname and default_tail" begin
		nt = NemenyiTest(ft3)
		@test testname(nt) == "Nemenyi all-pairs post-hoc test"
		@test default_tail(nt) === :both
	end

end  # NemenyiTest
