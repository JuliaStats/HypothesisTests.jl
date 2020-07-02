using HypothesisTests, Test

@testset "Clark-West tests" begin

    e1 = [ 9.78432,  12.73500,   8.67224,   2.62740,   5.60947,
           7.57809,   4.53774,   2.51084,   0.49290,   3.48394,
           9.46152,   9.41220,   2.36289,   0.34495,   3.33599,
           5.31357,   4.28219,  -3.74471,  -5.73575,  -3.71781 ]

    e2 =  [ 2.82053,   4.39754,  -1.78647,  -4.30662,   3.69526,
            3.37259,  -1.09002,  -0.50984,  -0.78973,   3.89077,
            7.52849,   2.82373,  -3.95838,  -0.13606,   4.54444,
            4.18216,   1.67993,  -5.38077,  -0.85686,   2.70479 ]

    atol = 1e-4
    cw_test = ClarkWestTest(e1, e2)
    @test pvalue(cw_test) ≈ 0.0002 atol=atol
    @test pvalue(cw_test, tail=:right) ≈ 0.0001 atol=atol

    cw_test = ClarkWestTest(e1, e2, 3)
    @test pvalue(cw_test) ≈ 0.0157 atol=atol
    @test pvalue(cw_test, tail=:right) ≈ 0.0079 atol=atol

    show(IOBuffer(), cw_test)

    @test_throws DimensionMismatch ClarkWestTest(rand(3), rand(4))

end