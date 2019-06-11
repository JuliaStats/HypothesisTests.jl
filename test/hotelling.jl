using HypothesisTests
using Test
using DelimitedFiles

@testset "Hotelling utility functions" begin
    HT = HypothesisTests

    @test HT.checkdims(rand(3, 2), rand(4, 2)) == (2, 3, 4)
    # Mismatched number of variables
    @test_throws DimensionMismatch HT.checkdims(rand(4, 3), rand(4, 6))
    # Empty observations
    @test_throws ArgumentError HT.checkdims(rand(0, 4), rand(2, 4))

    Sx = [ 3.0 -1.5 0.0
          -1.5  1.0 0.5
           0.0  0.5 1.0]
    Sy = fill(4.5, (3, 3))
    out = [3.5 0.5      1.5
           0.5 2.166667 1.833333
           1.5 1.833333 2.166667]
    @test HT.poolcov!(Sx, 2, Sy, 1) ≈ out atol=1e-6
    # Input is modified
    @test Sx ≈ out atol=1e-6

    # Positive semi-definite but not positive definite
    P = [ 1.0000  0.7426  0.1601 -0.7000 0.5500
          0.7426  1.0000 -0.2133 -0.5818 0.5000
          0.1601 -0.2133  1.0000 -0.1121 0.1000
         -0.7000 -0.5818 -0.1121  1.0000 0.4500
          0.5500  0.5000  0.1000  0.4500 1.0000]
    # Positive definite
    Q = [ 1.0 -0.5  0.0
         -0.5  1.0 -0.5
          0.0 -0.5  1.0]
    @test @inferred(HT.At_Binv_A(ones(5), P)) ≈ -0.8008792 atol=1e-6
    @test @inferred(HT.At_Binv_A(ones(3), Q)) ≈ 10.0
end

@testset "One sample Hotelling's T²" begin
    # Columns are calcium, iron, protein, vitamin A, vitamin C
    nutrient = readdlm(joinpath(@__DIR__, "data", "nutrient.txt"))[:,2:end]

    t = OneSampleHotellingT2Test(nutrient, [1000, 15, 60, 800, 75])
    @test nobs(t) == 737
    @test dof(t) == (5, 732)
    @test pvalue(t) ≈ 0.0 atol=eps()
    @test t.T² ≈ 1758.5413137 atol=1e-6
    @test t.F ≈ 349.7968048 atol=1e-6
    let out = sprint(show, t)
        @test occursin("reject h_0", out) && !occursin("fail to", out)
    end

    # Columns are survey answers: husbands' to questions 1-4, wives' to questions 1-4
    spouse = readdlm(joinpath(@__DIR__, "data", "spouse.txt"))

    # Paired
    p = OneSampleHotellingT2Test(spouse[:,1:4], spouse[:,5:end])
    @test nobs(p) == 30
    @test dof(p) == (4, 26)
    @test pvalue(p) ≈ 0.039369144 atol=1e-6
    @test p.T² ≈ 13.127840261 atol=1e-6
    @test p.F ≈ 2.942446955 atol=1e-6
    let out = sprint(show, p)
        @test occursin("reject h_0", out) && !occursin("fail to", out)
    end
end

@testset "Two sample Hotelling's T²" begin
    # Columns are type, length, left, right, bottom, top, diag
    swiss = readdlm(joinpath(@__DIR__, "data", "swiss3.txt"))
    genuine = convert(Matrix{Float64}, swiss[view(swiss, :, 1) .== "real", 2:end])
    counterfeit = convert(Matrix{Float64}, swiss[view(swiss, :, 1) .== "fake", 2:end])

    eq = EqualCovHotellingT2Test(genuine, counterfeit)
    @test nobs(eq) == (100, 100)
    @test dof(eq) == (6, 193)
    @test pvalue(eq) ≈ 0.0 atol=eps()
    @test eq.T² ≈ 2412.4506855 atol=1e-6
    @test eq.F ≈ 391.9217023 atol=1e-6
    @test occursin("reject h_0", sprint(show, eq))

    un = UnequalCovHotellingT2Test(genuine, counterfeit)
    @test nobs(un) == (100, 100)
    @test dof(un) == (6, 193)
    @test pvalue(un) ≈ 0.0 atol=eps()
    @test un.T² ≈ 2412.4506855 atol=1e-6
    @test un.F ≈ 391.9217023 atol=1e-6
    @test occursin("reject h_0", sprint(show, un))
end
