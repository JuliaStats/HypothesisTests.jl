module SimulateBinomial
    using Base.Test
    using Distributions, HypothesisTests

    include("utils.jl")

    has_ci(::Type{BinomialTest}) = true
    has_p_value(::Type{BinomialTest}) = true
    function call_test(::Type{BinomialTest}, x, simulation)
        return BinomialTest(sum(x), length(x), mean(simulation.d_x))
    end

    @testset "BinomialTest" begin
        results = simulate(
            OneSampleSimulation(Bernoulli(0.5), 100),
            BinomialTest,
            100_000,
        )
        analyze_simulation(results)
    end
end
