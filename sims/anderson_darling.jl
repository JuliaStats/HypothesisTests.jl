module SimulateAndersonDarling
    using Base.Test
    using Distributions, HypothesisTests

    include("utils.jl")

    has_ci(::Type{OneSampleADTest}) = false
    has_p_value(::Type{OneSampleADTest}) = true
    function call_test(::Type{OneSampleADTest}, x, simulation)
        return OneSampleADTest(x, simulation.d_x)
    end

    @testset "OneSampleADTest" begin
        results = simulate(
            OneSampleSimulation(Normal(0.0, 1.0), 100),
            OneSampleADTest,
            100_000,
        )
        analyze_simulation(results)
    end
end
