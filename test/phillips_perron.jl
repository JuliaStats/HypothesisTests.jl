using HypothesisTests, Test
using StableRNGs

@testset "Phillips-Perron" begin
    # A unit root process with MA(1) errors, and an AR(1) process around a linear trend
    rng = StableRNG(20261126)
    e = randn(rng, 201)
    y_unit = cumsum(e[2:end] .+ 0.5 .* e[1:end-1])
    u = randn(rng, 200)
    x = zeros(200)
    x[1] = u[1]
    for t in 2:200
        x[t] = 0.7 * x[t-1] + u[t]
    end
    y_stat = 10 .+ 0.05 .* (1:200) .+ x
    @test y_unit[1] ≈ -0.10706402692042377
    @test y_stat[1] ≈ 9.741406184825498

    # Reference values from Python arch 7.2.0,
    # PhillipsPerron(y, lags=lag, trend=..., test_type="tau"). arch 7.2.0 computes p-values
    # and critical values from the MacKinnon (1994, 2010) tables that ADFTest uses. arch
    # 9.0.0 gives the same statistics, but its p-values come from newer simulations.
    # R's urca 1.3.4 (ur.pp with type = "Z-tau") differs in the third decimal, because its
    # correction term uses the sample moments of y_t instead of the standard error of the
    # coefficient on y_{t-1}.
    cv_none = [-2.5770667644756444, -1.942430759336949, -1.615550641718986]
    cv_constant = [-3.4636447617687436, -2.8761761179270766, -2.57457158581854]
    cv_trend = [-4.004997848936356, -3.4327862452981046, -3.140144918368515]
    refs = [
        (y_unit, :none, 4, -0.44908578053904963, 0.5169254545425, cv_none),
        (y_unit, :constant, 4, -0.8257384350223911, 0.8114444489749435, cv_constant),
        (y_unit, :trend, 15, -1.4431932491667636, 0.8477835854068627, cv_trend),
        (y_stat, :none, 15, 1.1164153703669417, 0.9309576252537555, cv_none),
        (y_stat, :constant, 0, -2.5339882994395464, 0.10742206110458097, cv_constant),
        (y_stat, :trend, 4, -5.770384574812235, 6.581199014696416e-06, cv_trend),
    ]
    for (y, deterministic, lag, stat, p, cv) in refs
        t = PhillipsPerronTest(y; deterministic=deterministic, lag=lag)
        @test t.n == 199
        @test t.deterministic == deterministic
        @test t.lag == lag
        @test t.stat ≈ stat
        @test pvalue(t) ≈ p
        @test t.cv ≈ cv
    end

    @testset "Default lag" begin
        # floor(4 * (199 / 100)^(1 / 4)) = 4, as in urca::ur.pp(y, lags = "short")
        t = PhillipsPerronTest(y_unit)
        @test t.deterministic == :constant
        @test t.lag == 4
        @test t.stat ≈ -0.8257384350223911
        @test PhillipsPerronTest(y_stat; deterministic=:trend).lag == 4
        @test PhillipsPerronTest(randn(StableRNG(1), 50)).lag == 3
    end

    @testset "No lags gives the Dickey-Fuller test" begin
        for y in (y_unit, y_stat), deterministic in (:none, :constant, :trend)
            t = PhillipsPerronTest(y; deterministic=deterministic, lag=0)
            adf = ADFTest(y, deterministic, 0)
            @test t.stat ≈ adf.stat
            @test t.coef ≈ adf.coef
            @test pvalue(t) ≈ pvalue(adf)
        end
    end

    @testset "Invalid arguments" begin
        @test_throws ArgumentError PhillipsPerronTest(y_unit; deterministic=:squared_trend)
        @test_throws ArgumentError PhillipsPerronTest(y_unit; lag=-1)
        @test_throws ArgumentError PhillipsPerronTest(y_unit; lag=199)
        @test_throws ArgumentError PhillipsPerronTest([1.0, 2.0, 3.0]; deterministic=:trend)
    end

    @testset "Show" begin
        t = PhillipsPerronTest(y_unit)
        @test repr("text/plain", t; context=:compact => true) ==
            """
            Phillips-Perron unit root test
            ------------------------------
            Population details:
                parameter of interest:   coefficient on lagged non-differenced variable
                value under h_0:         0
                point estimate:          -0.00431964

            Test summary:
                outcome with 95% confidence: fail to reject h_0
                one-sided p-value:           0.8114

            Details:
                sample size in regression:          199
                lags in long-run variance:          4
                Z_t statistic:                      -0.825738
                Critical values at 1%, 5%, and 10%: -3.46364, -2.87618, -2.57457
            """
    end
end
