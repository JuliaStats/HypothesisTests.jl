using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "F-tests" begin

    @testset "Basic variance F-test" begin
        Random.seed!(12)
        y1_h0 = 4 .+ randn(500)
        y2_h0 = 4 .+ randn(400)

        t = VarianceFTest(y1_h0, y2_h0)

        @test t.n_x == 500
        @test t.n_y == 400
        @test t.df_x == 499
        @test t.df_y == 399
        @test t.F ≈ 0.974693 rtol = 1e-5
        @test pvalue(t) ≈ 0.784563 rtol = 1e-5
        @test default_tail(t) == :both

        y1_h1 = 0.7*randn(200)
        y2_h1 = 1.3*randn(120)

        t = VarianceFTest(y1_h1, y2_h1)

        @test t.n_x == 200
        @test t.n_y == 120
        @test t.df_x == 199
        @test t.df_y == 119
        @test t.F ≈ 0.367547 rtol = 1e-5
        @test pvalue(t) < 1e-8
        @test default_tail(t) == :both
        @test pvalue(t; tail = :left) < 1e-8
        @test pvalue(t; tail = :right) > 1.0 - 1e-8
    end

    @testset "Homoscedasticity" begin
        Random.seed!(12)
        y_h0 = randn(200)

        t = VarianceFTest(y_h0; firstlast = 100)

        @test t.n_x == t.n_y == 100
        @test t.df_x == t.df_y == 99
        @test t.F ≈ 1.191138 rtol = 1e-5
        @test pvalue(t) ≈ 0.385713 rtol = 1e-5
        @test default_tail(t) == :both

        y_h1 = [1.3*randn(100); 0.7*randn(100)]

        t = VarianceFTest(y_h1; firstlast = 80)

        @test t.n_x == t.n_y == 80
        @test t.df_x == t.df_y == 79
        @test t.F ≈ 3.965640 rtol = 1e-5
        @test pvalue(t) ≈ 3.972947e-9 rtol = 1e-5
        @test default_tail(t) == :both

        @test_throws ArgumentError VarianceFTest(y_h1; firstlast = 101)

        @test pvalue(t; tail = :left) > 1.0 - 1e-8
        @test pvalue(t; tail = :right) < 1e-8

        t = VarianceFTest(y_h0)
        @test t.n_x == t.n_y == 100
    end

end