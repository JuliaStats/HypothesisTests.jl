using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "Mann-Whitney" begin
    @testset "Basic exact test" begin
        test = ExactMannWhitneyUTest([1:10;], [2.1:2:21;])
        @test test.nx == 10
        @test test.ny == 10
        @test test.ranks ==
              [1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 3, 6, 9, 12, 15, 16, 17, 18, 19, 20]
        @test iszero(test.tie_adjustment)
        @test test.U == 20
        @test test.median == -5.6

        @test default_tail(test) == :both
        @test repr(test) == """
        Exact Mann-Whitney U test
        -------------------------
        Population details:
            parameter of interest:   Location parameter (pseudomedian)
            value under h_0:         0
            point estimate:          -5.6

        Test summary:
            outcome with 95% confidence: reject h_0
            two-sided p-value:           0.0232

        Details:
            number of observations in each group: [10, 10]
            Mann-Whitney-U statistic:             20.0
            rank sums:                            [75.0, 135.0]
            adjustment for ties:                  0.0
        """

        # Two-sided
        for kwargs in ((), (; tail=:both))
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2.1:2:21;]);
                                       kwargs...)) - 0.0232) <= 1e-4
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1:10;]);
                                       kwargs...)) - 0.0232) <= 1e-4
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1.5:10:100;], [2.1:2:21;]);
                                       kwargs...)) - 0.0068) <= 1e-4
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1.5:10:100;]);
                                       kwargs...)) - 0.0068) <= 1e-4
        end

        # Left tail
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2.1:2:21;]); tail=:left)) -
                  0.0116) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1:10;]); tail=:left)) -
                  0.9907) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1.5:10:100;], [2.1:2:21;]);
                                   tail=:left)) - 0.9974) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1.5:10:100;]);
                                   tail=:left)) - 0.0034) <= 1e-4

        # Right tail
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2.1:2:21;]);
                                   tail=:right)) - 0.9907) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1:10;]);
                                   tail=:right)) - 0.0116) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1.5:10:100;], [2.1:2:21;]);
                                   tail=:right)) - 0.0034) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1.5:10:100;]);
                                   tail=:right)) - 0.9974) <= 1e-4
    end

    @testset "Exact with ties" begin
        show(IOBuffer(), ExactMannWhitneyUTest([1:10;], [1:10;]))

        # Two-sided
        for kwargs in ((), (; tail=:both))
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:10;]); kwargs...)) -
                      1) <= 1e-4
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:11;]); kwargs...)) -
                      0.5096) <= 1e-4
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:11;], [1:10;]); kwargs...)) -
                      0.5096) <= 1e-4
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:5; ones(5)]);
                                       kwargs...)) - 0.0057) <= 1e-4
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:5; ones(5)], [1:10;]);
                                       kwargs...)) - 0.0057) <= 1e-4
        end

        # Left tail
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:10;]); tail=:left)) -
                  0.5296) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:11;]); tail=:left)) -
                  0.2548) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:11;], [1:10;]); tail=:left)) -
                  0.7634) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:5; ones(5)]);
                                   tail=:left)) - 0.9979) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:5; ones(5)], [1:10;]);
                                   tail=:left)) - 0.0028) <= 1e-4

        # Right tail
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:10;]); tail=:right)) -
                  0.5296) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:11;]); tail=:right)) -
                  0.7634) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:11;], [1:10;]); tail=:right)) -
                  0.2548) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:5; ones(5)]);
                                   tail=:right)) - 0.0028) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:5; ones(5)], [1:10;]);
                                   tail=:right)) - 0.9979) <= 1e-4
    end

    @testset "Exact with ties and unequal lengths" begin
        test1 = ExactMannWhitneyUTest([1:10;], [2:2:24;])
        test2 = ExactMannWhitneyUTest([2:2:24;], [1:10;])
        @test test1.nx == test2.ny == 10
        @test test2.nx == test1.ny == 12
        @test test1.ranks == test2.ranks[[13:22; 1:12]]
        @test test1.tie_adjustment == test2.tie_adjustment
        @test test1.U == 22.5
        @test test1.U + test2.U == 120
        @test test1.median == -test2.median

        @test repr(test1) == """
        Exact Mann-Whitney U test
        -------------------------
        Population details:
            parameter of interest:   Location parameter (pseudomedian)
            value under h_0:         0
            point estimate:          -7.5

        Test summary:
            outcome with 95% confidence: reject h_0
            two-sided p-value:           0.0120

        Details:
            number of observations in each group: [10, 12]
            Mann-Whitney-U statistic:             22.5
            rank sums:                            [77.5, 175.5]
            adjustment for ties:                  30.0
        """
        @test repr(test2) == """
        Exact Mann-Whitney U test
        -------------------------
        Population details:
            parameter of interest:   Location parameter (pseudomedian)
            value under h_0:         0
            point estimate:          7.5

        Test summary:
            outcome with 95% confidence: reject h_0
            two-sided p-value:           0.0120

        Details:
            number of observations in each group: [12, 10]
            Mann-Whitney-U statistic:             97.5
            rank sums:                            [175.5, 77.5]
            adjustment for ties:                  30.0
        """

        # Two-sided
        for kwargs in ((), (; tail=:both))
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:2:24;]);
                                       kwargs...)) - 0.0120) <= 1e-4
            @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:2:24;], [1:10;]);
                                       kwargs...)) - 0.0120) <= 1e-4
        end

        # Left tail
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:2:24;]); tail=:left)) -
                  0.0060) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:2:24;], [1:10;]); tail=:left)) -
                  0.9949) <= 1e-4

        # Right tail
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:2:24;]); tail=:right)) -
                  0.9949) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:2:24;], [1:10;]); tail=:right)) -
                  0.0060) <= 1e-4
    end

    @testset "Approximate test" begin
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;]))) -
                  0.0257) <= 1e-4
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([2.1:2:21;], [1:10;]))) -
                  0.0257) <= 1e-4
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1.5:10:100;], [2.1:2:21;]))) -
                  0.0091) <= 1e-4
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([2.1:2:21;], [1.5:10:100;]))) -
                  0.0091) <= 1e-4
        @test default_tail(ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;])) == :both
        show(IOBuffer(), ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;]))
    end

    @testset "Approximate with ties" begin
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:10;], [1:10;]))) - 1) <=
              1e-4
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:10;], [2:11;]))) - 0.4948) <=
              1e-4
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([2:11;], [1:10;]))) - 0.4948) <=
              1e-4
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:10;], [1:5; ones(5)]))) -
                  0.0076) <= 1e-4
        @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:5; ones(5)], [1:10;]))) -
                  0.0076) <= 1e-4
        show(IOBuffer(), ApproximateMannWhitneyUTest([1:10;], [1:10;]))
    end

    @testset "Tests for automatic selection" begin
        @test abs(@inferred(pvalue(MannWhitneyUTest([1:10;], [2.1:2:21;]))) - 0.0232) <=
              1e-4
        @test abs(@inferred(pvalue(MannWhitneyUTest([1:10;], [2:11;]))) - 0.4948) <= 1e-4
        show(IOBuffer(), MannWhitneyUTest([1:10;], [2.1:2:21;]))
    end

    @testset "Issue #39" begin
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest(Float32[1:10;], Float32[2:11;]))) -
                  0.5096) <= 1e-4
    end

    @testset "Issue #113" begin
        @test abs(pvalue(ApproximateMannWhitneyUTest(Float32[1:10;], Float32[2:11;])) -
                  0.4948) <= 1e-4
    end

    @testset "Issue #126 pvalues above 1" begin
        A = [1.34937, 1.75722, 0.276514, 1.04546, 1.69085, 0.738085, 2.36313]
        B = [2.62325, 1.16533, 1.1327, 0.728714]

        m = ExactMannWhitneyUTest(A, B)
        p = @inferred(pvalue(m; tail=:both))
        @test p == 1

        A = [12, 10, 7, 6, 3, 1]
        B = [11, 9, 8, 5, 4, 2]

        m = MannWhitneyUTest(A, B)
        p = @inferred(pvalue(m; tail=:both))

        @test p == 1

        m = ApproximateMannWhitneyUTest(A, B)
        p = @inferred(pvalue(m; tail=:both))
    end
end
