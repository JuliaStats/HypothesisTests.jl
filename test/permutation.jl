using HypothesisTests, Test, Random, Statistics

@testset "Permutation" begin
@test isapprox(pvalue(PermutationTest(mean,[1,2,3],[4,5,6]), tail=:both), 0.1)
@test isapprox(pvalue(PermutationTest(mean,[4,5,6],[1,2,3]), tail=:both), 0.1)
@test isapprox(pvalue(PermutationTest(mean,[1,2,3],[4,5,6]), tail=:left), 0.05)
@test isapprox(pvalue(PermutationTest(mean,[4,5,6],[1,2,3]), tail=:right), 0.05)

Random.seed!(12345)
@test isapprox(pvalue(PermutationTest(mean,[1,2,3],[4,5,6],100), tail=:both), 0.09)
@test isapprox(pvalue(PermutationTest(mean,[4,5,6],[1,2,3],100), tail=:both), 0.04)
@test isapprox(pvalue(PermutationTest(mean,[1,2,3],[4,5,6],100), tail=:left), 0.08)
@test isapprox(pvalue(PermutationTest(mean,[4,5,6],[1,2,3],100), tail=:right), 0.07)
end
