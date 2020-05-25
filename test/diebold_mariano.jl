using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "Diebold-Mariano tests" begin
    e_ets = [0.002241270, -0.000661906, -1.200002613, -4.000218178, 3.799373533,
             3.401005111, -1.199692125, -0.800490048, -1.200048999, 3.399847229,
             7.000708948,  2.200989241, -4.800165569, -1.201043375, 3.400167556, 
             3.000709085,  0.400267950, -6.800169207, -2.401259050, 1.200112198]

    e_arima = [0.1209998,  8.0622269, -0.9189429, -1.7952436, 5.3431277, 2.0459823,
               -0.6073261, 0.5516816, -0.7633238, 4.1229472, 5.8785333, 1.2828331,
               -3.4854233, 1.1971370, 3.1709314, 2.3408457, 0.4962286, -5.9841362,
                0.1537114, 0.4207951]

    atol = 1e-3
    dm_test = DieboldMarianoTest(e_ets, e_arima)
    @test pvalue(dm_test) ≈ 0.8871 atol=atol
    @test pvalue(dm_test, tail=:right) ≈ 0.5565 atol=atol
    @test pvalue(dm_test, tail=:left) ≈ 0.4435 atol=atol

    dm_test = DieboldMarianoTest(e_ets, e_arima; lookahead=10)
    @test pvalue(dm_test) ≈ 0.9362 atol=atol
    @test pvalue(dm_test, tail=:right) ≈ 0.5319 atol=atol
    @test pvalue(dm_test, tail=:left) ≈ 0.4681 atol=atol

    dm_test = DieboldMarianoTest(e_ets, e_arima; loss=abs)
    @test pvalue(dm_test) ≈ 0.7818 atol=atol
    @test pvalue(dm_test, tail=:right) ≈ 0.3909 atol=atol
    @test pvalue(dm_test, tail=:left) ≈ 0.6091 atol=atol
    show(IOBuffer(), dm_test)

    @test_throws DimensionMismatch DieboldMarianoTest(rand(3), rand(4))
end
