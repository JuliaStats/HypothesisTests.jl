using HypothesisTests, Test

mutable struct TestTest <: HypothesisTests.HypothesisTest end

@testset "Common" begin
@test_throws DimensionMismatch HypothesisTests.check_same_length([1], [])
@test_throws ArgumentError HypothesisTests.check_level(1.0)

result = HypothesisTests.population_param_of_interest(TestTest())
@test result[1] == "not implemented yet"
@test isnan(result[2])
end
