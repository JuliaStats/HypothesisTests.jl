using HypothesisTests
using Test
using StatsBase
using DelimitedFiles

@testset "Equality of Variances" begin
    @testset "One-way ANOVA" begin
        # https://en.wikipedia.org/wiki/One-way_analysis_of_variance#Example
        groups = [
            [6, 8, 4, 5, 3, 4],
            [8, 12, 9, 11, 6, 8],
            [13, 9, 11, 8, 7, 12]
        ]
        t = OneWayANOVATest(groups)
        @test nobs(t) == fill(6, 3)
        @test dof(t) == (2,15)
        @test pvalue(t) ≈ 0.002 atol=1e-3
        @test occursin("reject h_0", sprint(show, t))

        # test splatting version
        t2 = OneWayANOVATest(groups...)
        @test nobs(t2) == nobs(t)
        @test dof(t2) == dof(t)
        @test pvalue(t2) == pvalue(t)
        @test HypothesisTests.teststatistic(t2) == HypothesisTests.teststatistic(t)
        @test sprint(show, t2) == sprint(show, t)

        show(IOContext(IOBuffer(), :table => true), t)
        show(IOBuffer(), t)

        # http://www.real-statistics.com/one-way-analysis-of-variance-anova/confidence-interval-anova/
        groups = [
            [51, 87, 50, 48, 79, 61, 53],
            [82, 91, 92, 80, 52, 79, 73, 74],
            [79, 84, 74, 98, 63, 83, 85, 58],
            [85, 80, 65, 71, 67, 51],
        ]
        t = OneWayANOVATest(groups)
        @test nobs(t) == [7, 8, 8, 6]
        @test dof(t) == (3, 25)
        @test pvalue(t) ≈ 0.07276 atol=1e-6
        @test HypothesisTests.teststatistic(t) ≈ 2.62311 atol=1e-6
        @test occursin("reject h_0", sprint(show, t))
    end

    # http://www.real-statistics.com/one-way-analysis-of-variance-anova/homogeneity-variances/levenes-test/
    groups = [
        [51, 87, 50, 48, 79, 61, 53, 54],
        [82, 91, 92, 80, 52, 79, 73, 74],
        [79, 84, 74, 98, 63, 83, 85, 58],
        [85, 80, 65, 71, 67, 51, 63, 93],
    ]
    @testset "Levene" begin
        # with means
        l = LeveneTest(groups; statistic=mean)
        @test nobs(l) == fill(8, 4)
        @test dof(l) == (3,28)
        @test pvalue(l) ≈ 0.90357 atol=1e-4

        l2 = LeveneTest(groups...; statistic=mean)
        @test nobs(l2) == nobs(l)
        @test dof(l2) == dof(l)
        @test pvalue(l2) == pvalue(l)
        @test HypothesisTests.teststatistic(l2) == HypothesisTests.teststatistic(l)
        @test sprint(show, l2) == sprint(show, l)

        # with medians
        l = LeveneTest(groups; statistic=median)
        @test pvalue(l) ≈ 0.97971 atol=1e-4
        # with 10% trimmed means
        l = LeveneTest(groups; statistic=v -> mean(trim(v, prop=0.1)))
        @test pvalue(l) ≈ 0.90357 atol=1e-4
    end

    @testset "Fligner-Killeen" begin
        t = FlignerKilleenTest(groups)
        @test nobs(t) == fill(8, 4)
        @test dof(t) == 3
        @test pvalue(t) ≈ 0.9878 atol=1e-4
        @test HypothesisTests.teststatistic(t) ≈ 0.1311 atol=1e-5
        @test occursin("reject h_0", sprint(show, t))

        # test splatting version
        t2 = FlignerKilleenTest(groups...)
        @test nobs(t2) == nobs(t)
        @test dof(t2) == dof(t)
        @test pvalue(t2) == pvalue(t)
        @test HypothesisTests.teststatistic(t2) == HypothesisTests.teststatistic(t)
        @test sprint(show, t2) == sprint(show, t)
    end

    @testset "Brown-Forsythe" begin
        # https://www.itl.nist.gov/div898/handbook/eda/section3/eda35a.htm
        # Columns are gear diameter and batch number
        gear = readdlm(joinpath(@__DIR__, "data", "gear.txt"))
        samples = reshape(gear[:, 1], :, 10)

        l = BrownForsytheTest(collect(eachcol(samples)))
        @test nobs(l) == fill(10, 10)
        @test dof(l) == (9, 90)
        @test HypothesisTests.teststatistic(l) ≈ 1.705910 atol=1e-5
        @test pvalue(l) ≈ 0.0991 atol=1e-4
        @test occursin("reject h_0", sprint(show, l))

        l2 = BrownForsytheTest(eachcol(samples)...)
        @test nobs(l2) == nobs(l)
        @test dof(l2) == dof(l)
        @test HypothesisTests.teststatistic(l2) == HypothesisTests.teststatistic(l)
        @test pvalue(l2) == pvalue(l)
        @test sprint(show, l2) == sprint(show, l)
    end
end
