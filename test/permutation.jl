using HypothesisTests, Test

@testset "Permutation" begin
@test isapprox(pvalue(ExactPermutationTest([1,2,3],[4,5,6],mean), tail=:both), 0.1)
@test isapprox(pvalue(ExactPermutationTest([4,5,6],[1,2,3],mean), tail=:both), 0.1)
@test isapprox(pvalue(ExactPermutationTest([1,2,3],[4,5,6],mean), tail=:left), 0.05)
@test isapprox(pvalue(ExactPermutationTest([4,5,6],[1,2,3],mean), tail=:right), 0.05)

Random.seed!(12345)
@test isapprox(pvalue(ApproximatePermutationTest([1,2,3],[4,5,6],mean,100), tail=:both), 0.09)
@test isapprox(pvalue(ApproximatePermutationTest([4,5,6],[1,2,3],mean,100), tail=:both), 0.04)
@test isapprox(pvalue(ApproximatePermutationTest([1,2,3],[4,5,6],mean,100), tail=:left), 0.08)
@test isapprox(pvalue(ApproximatePermutationTest([4,5,6],[1,2,3],mean,100), tail=:right), 0.07)
end
