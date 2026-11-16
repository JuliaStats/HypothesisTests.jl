module SimulateMannWhitney
    using Base.Test
    using Distributions, HypothesisTests

    include("utils.jl")

    has_ci(::Type{ExactMannWhitneyUTest}) = false
    has_p_value(::Type{ExactMannWhitneyUTest}) = true
    function call_test(::Type{ExactMannWhitneyUTest}, x, y, simulation)
        return ExactMannWhitneyUTest(x, y)
    end

    has_ci(::Type{ApproximateMannWhitneyUTest}) = false
    has_p_value(::Type{ApproximateMannWhitneyUTest}) = true
    function call_test(::Type{ApproximateMannWhitneyUTest}, x, y, simulation)
        return ApproximateMannWhitneyUTest(x, y)
    end

    @testset "ExactMannWhitneyUTest" begin
        results = simulate(
            TwoSampleSimulation(Normal(0.0, 1.0), Normal(0.0, 1.0), 50, 100),
            ExactMannWhitneyUTest,
            100_000,
        )
        analyze_simulation(results)
    end

    @testset "ApproximateMannWhitneyUTest" begin
        results = simulate(
            TwoSampleSimulation(Normal(0.0, 1.0), Normal(0.0, 1.0), 100, 100),
            ApproximateMannWhitneyUTest,
            100_000,
        )
        analyze_simulation(results)
    end
end
