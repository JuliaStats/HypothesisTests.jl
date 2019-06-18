using HypothesisTests
using Test
using StatsBase
using DelimitedFiles

@testset "Bartlett's test" begin
    # Columns are type, length, left, right, bottom, top, diag
    swiss = readdlm(joinpath(@__DIR__, "data", "swiss3.txt"))
    genuine = convert(Matrix{Float64}, swiss[view(swiss, :, 1) .== "real", 2:end])
    counterfeit = convert(Matrix{Float64}, swiss[view(swiss, :, 1) .== "fake", 2:end])

    b = BartlettTest(genuine, counterfeit)
    @test nobs(b) == (100, 100)
    @test dof(b) == 21
    @test pvalue(b) ≈ 0.0 atol=1e-10
    @test b.L′ ≈ 121.8991235 atol=1e-6
    @test occursin("reject h_0", sprint(show, b))
end
