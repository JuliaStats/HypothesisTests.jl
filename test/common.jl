using HypothesisTests, Test

mutable struct TestTest <: HypothesisTests.HypothesisTest end

struct RefusingCITest <: HypothesisTests.HypothesisTest end
struct WorkingCITest <: HypothesisTests.HypothesisTest end
struct BrokenCITest <: HypothesisTests.HypothesisTest end
struct ArgumentBugCITest <: HypothesisTests.HypothesisTest end
const StubCITest = Union{RefusingCITest,WorkingCITest,BrokenCITest,ArgumentBugCITest}
HypothesisTests.confint(::RefusingCITest; level::Real=0.95) =
    throw(HypothesisTests.ComputationTooLarge("refusing to form too many pairwise estimates"))
HypothesisTests.confint(::WorkingCITest; level::Real=0.95) = (-1.0, 1.0)
HypothesisTests.confint(::BrokenCITest; level::Real=0.95) = error("a bug, not a refusal")
# The case that motivated giving refusals a type of their own. `confint` is generic over
# `HypothesisTest`, so while the refusal was an `ArgumentError` this method's bug was
# indistinguishable from one and was printed away instead of raising.
HypothesisTests.confint(::ArgumentBugCITest; level::Real=0.95) =
    throw(ArgumentError("a bug in someone else's confint, not a refusal"))
HypothesisTests.pvalue(::StubCITest; tail=:both) = 0.5
HypothesisTests.testname(::StubCITest) = "Stub"

@testset "Common" begin
@test_throws DimensionMismatch HypothesisTests.check_same_length([1], [])
@test_throws ArgumentError HypothesisTests.check_level(1.0)

result = HypothesisTests.population_param_of_interest(TestTest())
@test result[1] == "not implemented yet"
@test isnan(result[2])

# The interval `show` prints is a convenience of the display. A `confint` that refuses
# the sample it was handed costs that one line; any other failure is a bug and propagates.
@test HypothesisTests.show_confint(TestTest()) === nothing
@test HypothesisTests.show_confint(RefusingCITest()) === nothing
@test HypothesisTests.show_confint(WorkingCITest()) == (-1.0, 1.0)
@test_throws ErrorException HypothesisTests.show_confint(BrokenCITest())
@test_throws ArgumentError HypothesisTests.show_confint(ArgumentBugCITest())
@test !occursin("confidence interval", sprint(show, RefusingCITest()))
@test occursin("95% confidence interval: (-1.0, 1.0)", sprint(show, WorkingCITest()))
@test_throws ArgumentError sprint(show, ArgumentBugCITest())
end
