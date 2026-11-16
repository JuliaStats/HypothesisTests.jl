immutable OneSampleSimulation{T}
    d_x::T
    n_x::Int
end

immutable TwoSampleSimulation{S, T}
    d_x::S
    d_y::T
    n_x::Int
    n_y::Int
end

immutable SimulationAnalysis
    coverage::Vector{Bool}
    p_values::Vector{Float64}
    has_ci::Bool
    has_p_value::Bool
end

covers(ci, truth) = ci[1] <= truth <= ci[2]

function initialize(simulation::OneSampleSimulation)
    return rand(simulation.d_x, simulation.n_x)
end

function initialize(simulation::TwoSampleSimulation)
    return (
        rand(simulation.d_x, simulation.n_x),
        rand(simulation.d_y, simulation.n_y),
    )
end

function fill!(x, simulation::OneSampleSimulation)
    rand!(simulation.d_x, x)
    return
end

function fill!(x, y, simulation::TwoSampleSimulation)
    rand!(simulation.d_x, x)
    rand!(simulation.d_y, y)
    return
end

function simulate(simulation::OneSampleSimulation, test, n_sims)
    # Initialize core results.
    coverage = Array(Bool, n_sims)
    p_values = Array(Float64, n_sims)

    # Initialize sample.
    x = initialize(simulation)

    # Run simulations.
    for sim in 1:n_sims
        fill!(x, simulation)
        test_results = call_test(test, x, simulation)
        if has_ci(test)
            coverage[sim] = covers(ci(test_results), mean(simulation.d_x))
        end
        if has_p_value(test)
            p_values[sim] = pvalue(test_results)
        end
    end

    return SimulationAnalysis(
        coverage,
        p_values,
        has_ci(test),
        has_p_value(test),
    )
end

function simulate(simulation::TwoSampleSimulation, test, n_sims)
    # Initialize core results.
    coverage = Array(Bool, n_sims)
    p_values = Array(Float64, n_sims)

    # Initialize sample.
    x, y = initialize(simulation)

    # Run simulations.
    for sim in 1:n_sims
        fill!(x, y, simulation)
        test_results = call_test(test, x, y, simulation)
        if has_ci(test)
            coverage[sim] = covers(
                ci(test_results),
                mean(simulation.d_x) - mean(simulation.d_y),
            )
        end
        if has_p_value(test)
            p_values[sim] = pvalue(test_results)
        end
    end

    return SimulationAnalysis(
        coverage,
        p_values,
        has_ci(test),
        has_p_value(test),
    )
end

# τ: Inverted tolerance for failure. Set to 12, which implies a 12-sigma test,
# which should fail during 1.776482112077648e-33 simulation tests. Even when
# we run many tests, this level of failure is rare -- moreover, the SEM is
# also small for large n, so the tests are quite stringent when the hypothesis
# is false.
function analyze_simulation(results::SimulationAnalysis, τ = 12)
    if results.has_ci
        n_sims = length(results.coverage)
        sem_coverage = 0.95 * (1 - 0.95) / sqrt(n_sims)
        @test abs(mean(results.coverage) - 0.95) <= τ * sem_coverage
    end

    if results.has_p_value
        n_sims = length(results.p_values)

        sem_mean = sqrt(var(Uniform(0, 1)) / n_sims)
        @test abs(mean(results.p_values) - 0.50) <= τ * sem_mean

        for α in (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
            sem_quantile = sqrt(var(Bernoulli(α)) / n_sims)
            @test abs(mean(results.p_values .<= α) - α) <= τ * sem_quantile
        end
    end

    return
end
