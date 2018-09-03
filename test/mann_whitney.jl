using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "Mann-Whitney" begin
@testset "Basic exact test" begin
	@test abs(pvalue(ExactMannWhitneyUTest([1:10;], [2.1:2:21;])) - 0.0232) <= 1e-4
	@test abs(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1:10;])) - 0.0232) <= 1e-4
	@test abs(pvalue(ExactMannWhitneyUTest([1.5:10:100;], [2.1:2:21;])) - 0.0068) <= 1e-4
	@test abs(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1.5:10:100;])) - 0.0068) <= 1e-4
	@test default_tail(ExactMannWhitneyUTest([1:10;], [2.1:2:21;])) == :both
	show(IOBuffer(), ExactMannWhitneyUTest([1:10;], [2.1:2:21;]))
end

@testset "Exact with ties" begin
	@test abs(pvalue(ExactMannWhitneyUTest([1:10;], [1:10;])) - 1) <= 1e-4
	@test abs(pvalue(ExactMannWhitneyUTest([1:10;], [2:11;])) - 0.5096) <= 1e-4
	@test abs(pvalue(ExactMannWhitneyUTest([2:11;], [1:10;])) - 0.5096) <= 1e-4
	@test abs(pvalue(ExactMannWhitneyUTest([1:10;], [1:5; ones(5)])) - 0.0057) <= 1e-4
	@test abs(pvalue(ExactMannWhitneyUTest([1:5; ones(5)], [1:10;])) - 0.0057) <= 1e-4
	show(IOBuffer(), ExactMannWhitneyUTest([1:10;], [1:10;]))
end

@testset "Exact with ties and unequal lengths" begin
	@test abs(pvalue(ExactMannWhitneyUTest([1:10;], [2:2:24;])) - 0.0118) <= 1e-4
	@test abs(pvalue(ExactMannWhitneyUTest([2:2:24;], [1:10;])) - 0.0118) <= 1e-4
	show(IOBuffer(), ExactMannWhitneyUTest([1:10;], [2:2:24;]))
end

@testset "Approximate test" begin
	@test abs(pvalue(ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;])) - 0.0257) <= 1e-4
	@test abs(pvalue(ApproximateMannWhitneyUTest([2.1:2:21;], [1:10;])) - 0.0257) <= 1e-4
	@test abs(pvalue(ApproximateMannWhitneyUTest([1.5:10:100;], [2.1:2:21;])) - 0.0091) <= 1e-4
	@test abs(pvalue(ApproximateMannWhitneyUTest([2.1:2:21;], [1.5:10:100;])) - 0.0091) <= 1e-4
	@test default_tail(ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;])) == :both
	show(IOBuffer(), ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;]))
end

@testset "Approximate with ties" begin
	@test abs(pvalue(ApproximateMannWhitneyUTest([1:10;], [1:10;])) - 1) <= 1e-4
	@test abs(pvalue(ApproximateMannWhitneyUTest([1:10;], [2:11;])) - 0.4948) <= 1e-4
	@test abs(pvalue(ApproximateMannWhitneyUTest([2:11;], [1:10;])) - 0.4948) <= 1e-4
	@test abs(pvalue(ApproximateMannWhitneyUTest([1:10;], [1:5; ones(5)])) - 0.0076) <= 1e-4
	@test abs(pvalue(ApproximateMannWhitneyUTest([1:5; ones(5)], [1:10;])) - 0.0076) <= 1e-4
	show(IOBuffer(), ApproximateMannWhitneyUTest([1:10;], [1:10;]))
end

@testset "Tests for automatic selection" begin
	@test abs(pvalue(MannWhitneyUTest([1:10;], [2.1:2:21;])) - 0.0232) <= 1e-4
	@test abs(pvalue(MannWhitneyUTest([1:10;], [2:11;])) - 0.4948) <= 1e-4
	show(IOBuffer(), MannWhitneyUTest([1:10;], [2.1:2:21;]))
end

@testset "Issue #39" begin
	@test abs(pvalue(ExactMannWhitneyUTest(Float32[1:10;], Float32[2:11;])) - 0.5096) <= 1e-4
end
end
