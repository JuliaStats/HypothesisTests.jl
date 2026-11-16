module SimulateTTest
    using Base.Test
    using Distributions, HypothesisTests

    include("utils.jl")

    has_ci(::Type{OneSampleTTest}) = true
    has_p_value(::Type{OneSampleTTest}) = true
    function call_test(::Type{OneSampleTTest}, x, simulation)
        return OneSampleTTest(x)
    end

    has_ci(::Type{EqualVarianceTTest}) = true
    has_p_value(::Type{EqualVarianceTTest}) = true
    function call_test(::Type{EqualVarianceTTest}, x, y, simulation)
        return EqualVarianceTTest(x, y)
    end

    has_ci(::Type{UnequalVarianceTTest}) = true
    has_p_value(::Type{UnequalVarianceTTest}) = true
    function call_test(::Type{UnequalVarianceTTest}, x, y, simulation)
        return UnequalVarianceTTest(x, y)
    end

    @testset "OneSampleTTest" begin
        results = simulate(
            OneSampleSimulation(Normal(0.0, 1.0), 100),
            OneSampleTTest,
            100_000,
        )
        analyze_simulation(results)
    end

    @testset "EqualVarianceTTest" begin
        results = simulate(
            TwoSampleSimulation(Normal(0.0, 1.0), Normal(0.0, 1.0), 100, 100),
            EqualVarianceTTest,
            100_000,
        )
        analyze_simulation(results)
    end

    @testset "UnequalVarianceTTest" begin
        results = simulate(
            TwoSampleSimulation(Normal(0.0, 1.0), Normal(0.0, 2.0), 100, 100),
            UnequalVarianceTTest,
            100_000,
        )
        analyze_simulation(results)
    end
end
