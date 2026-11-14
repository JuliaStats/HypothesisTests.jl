module SimulateZTest
    # TODO: We get occasional failures here. We should consider forcing
    # OneSampleZTest to provide the true variance rather than allow it to
    # estimate it and pretend the estimate is known with certainty. You
    # reliably detect these failtures if you simulate 1,000,000 times instead
    # of 10,000. See t.jl for a contrast where all assumptions are satisfied.

    using Base.Test
    using Distributions, HypothesisTests

    include("utils.jl")

    has_ci(::Type{OneSampleZTest}) = true
    has_p_value(::Type{OneSampleZTest}) = true
    function call_test(::Type{OneSampleZTest}, x, simulation)
        return OneSampleZTest(x)
    end

    has_ci(::Type{EqualVarianceZTest}) = true
    has_p_value(::Type{EqualVarianceZTest}) = true
    function call_test(::Type{EqualVarianceZTest}, x, y, simulation)
        return EqualVarianceZTest(x, y)
    end

    has_ci(::Type{UnequalVarianceZTest}) = true
    has_p_value(::Type{UnequalVarianceZTest}) = true
    function call_test(::Type{UnequalVarianceZTest}, x, y, simulation)
        return UnequalVarianceZTest(x, y)
    end

    @testset "OneSampleZTest" begin
        results = simulate(
            OneSampleSimulation(Normal(0.0, 1.0), 100),
            OneSampleZTest,
            10_000,
        )
        analyze_simulation(results)
    end

    @testset "EqualVarianceZTest" begin
        results = simulate(
            TwoSampleSimulation(Normal(0.0, 1.0), Normal(0.0, 1.0), 100, 100),
            EqualVarianceZTest,
            10_000,
        )
        analyze_simulation(results)
    end

    @testset "UnequalVarianceZTest" begin
        results = simulate(
            TwoSampleSimulation(Normal(0.0, 1.0), Normal(0.0, 2.0), 100, 100),
            UnequalVarianceZTest,
            10_000,
        )
        analyze_simulation(results)
    end
end
