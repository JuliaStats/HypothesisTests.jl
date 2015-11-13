using HypothesisTests, Base.Test

@test_throws DimensionMismatch HypothesisTests.check_same_length([1], [])
@test_throws ArgumentError HypothesisTests.check_alpha(0.0)

type TestTest <: HypothesisTests.HypothesisTest end
result = HypothesisTests.population_param_of_interest(TestTest())
@test result[1] == "not implemented yet"
@test isnan(result[2])
