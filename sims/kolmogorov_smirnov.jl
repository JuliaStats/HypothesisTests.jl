module SimulateKolmogorovSmirnov
    using Base.Test
    using Distributions, HypothesisTests

    include("utils.jl")

    has_ci(::Type{ExactOneSampleKSTest}) = false
    has_p_value(::Type{ExactOneSampleKSTest}) = true
    function call_test(::Type{ExactOneSampleKSTest}, x, simulation)
        return ExactOneSampleKSTest(x, simulation.d_x)
    end

    @testset "ExactOneSampleKSTest" begin
        results = simulate(
            OneSampleSimulation(Normal(0.0, 1.0), 100),
            ExactOneSampleKSTest,
            100_000,
        )
        analyze_simulation(results)
    end
end
