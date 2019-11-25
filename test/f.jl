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
        @test t.F ≈ 0.974693 rtol=1e-5
        @test pvalue(t) ≈ 0.784563 rtol=1e-5
        @test pvalue(t, tail=:left) ≈ 0.392281 rtol=1e-5
        @test pvalue(t, tail=:right) ≈ 0.607718 rtol=1e-5
        @test default_tail(t) == :both

        t = VarianceFTest(y2_h0, y1_h0)

        @test t.n_x == 400
        @test t.n_y == 500
        @test t.df_x == 399
        @test t.df_y == 499
        @test t.F ≈ 1.025963 rtol=1e-5
        @test pvalue(t) ≈ 0.784563 rtol=1e-5
        @test pvalue(t, tail=:right) ≈ 0.392281 rtol=1e-5
        @test pvalue(t, tail=:left) ≈ 0.607718 rtol=1e-5
        @test default_tail(t) == :both

        y1_h1 = 0.7*randn(200)
        y2_h1 = 1.3*randn(120)

        t = VarianceFTest(y1_h1, y2_h1)

        @test t.n_x == 200
        @test t.n_y == 120
        @test t.df_x == 199
        @test t.df_y == 119
        @test t.F ≈ 0.367547 rtol=1e-5
        @test pvalue(t) < 1e-8
        @test default_tail(t) == :both
        @test pvalue(t, tail=:left) < 1e-8
        @test pvalue(t, tail=:right) > 1.0 - 1e-8
    end
end
