module SimulateKruskalWallis
    using Base.Test
    using Distributions, HypothesisTests

    include("utils.jl")

    has_ci(::Type{KruskalWallisTest}) = false
    has_p_value(::Type{KruskalWallisTest}) = true
    function call_test(::Type{KruskalWallisTest}, x, y, simulation)
        return KruskalWallisTest(x, y)
    end

    @testset "KruskalWallisTest" begin
        results = simulate(
            TwoSampleSimulation(Normal(0.0, 1.0), Normal(0.0, 1.0), 100, 100),
            KruskalWallisTest,
            100_000,
        )
        analyze_simulation(results)
    end
end
