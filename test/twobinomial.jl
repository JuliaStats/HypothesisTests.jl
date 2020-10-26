using HypothesisTests, Test

@testset "Two Binomial" begin

    #https://www.lexjansen.com/wuss/2016/127_Final_Paper_PDF.pdf
    #Constructing Confidence Intervals for the Differences of Binomial Proportions in SAS® Will Garner, Gilead Sciences, Inc., Foster City, CA
    t = TwoSampleBinomialTest(56, 70, 48, 80, -0.1; ptype = :difference, htype = :superiority, alpha = 0.025,  method = :mn)
    @test_ci_approx t.ci (0.052829713256921054, 0.33817294018078103)

    t = TwoSampleBinomialTest(56, 70, 48, 80, -0.1; ptype = :difference, htype = :superiority, alpha = 0.025,  method=:mnmee)
    @test_ci_approx t.ci (0.053333849176126207, 0.33772946376708673)

    t = TwoSampleBinomialTest(56, 70, 48, 80, -0.1; ptype = :difference, htype = :superiority, alpha = 0.025,  method=:nhscc)
    @test_ci_approx t.ci (0.04276787200288948, 0.3421862784381806)

    #Recommended confidence intervals for two independent binomial proportions DOI: 10.1177/0962280211415469
    t = TwoSampleBinomialTest(7, 34, 1, 34, -0.1; ptype = :difference, htype =:superiority, alpha=0.025,  method=:nhs)
    @test_ci_approx t.ci (0.018921443885772632, 0.3403686870327077)

    t = TwoSampleBinomialTest(56, 70, 48, 80, -0.1; ptype = :difference, htype =:superiority, alpha=0.025,  method=:wald)
    @test_ci_approx t.ci (0.05750489914903689, 0.34249510085096324)

    t = TwoSampleBinomialTest(56, 70, 48, 80, -0.1; ptype = :difference, htype =:superiority, alpha=0.025,  method=:waldcc)
    @test_ci_approx t.ci (0.04411204200617975, 0.3558879579938204)

    t = TwoSampleBinomialTest(56, 70, 48, 80, -0.1; ptype = :difference, htype =:superiority, alpha=0.025,  method=:ac)
    @test_ci_approx t.ci (0.052452926638053426, 0.33575845547576766)

    t = TwoSampleBinomialTest(56, 70, 48, 80, -0.1; ptype = :difference, htype =:superiority, alpha=0.025,  method=:ha)
    @test_ci_approx t.ci (0.042263997307960856, 0.35773600269203926)

    t = TwoSampleBinomialTest(56, 70, 48, 80, -0.1; ptype = :difference, htype =:superiority, alpha=0.025,  method=:mover)
    @test_ci_approx t.ci (0.05243147240236476, 0.33387265403690625)

    #https://rdrr.io/cran/ORCI/man/Woolf.CI.html
    t = TwoSampleBinomialTest(2, 14, 1, 11, 0.1; ptype = :oddratio, htype =:superiority, alpha=0.025,  method=:woolf)
    @test_ci_approx t.ci (0.13106039368845537, 21.194639353677303)

    t = TwoSampleBinomialTest(2, 14, 1, 11, 0.1; ptype = :oddratio, htype =:superiority, alpha=0.025,  method=:mn)
    @test_ci_approx t.ci (0.17631575652726034, 14.876663936201446)

    t = TwoSampleBinomialTest(2, 14, 1, 11, 0.1; ptype = :oddratio, htype =:superiority, alpha=0.025,  method=:awoolf)
    @test_ci_approx t.ci (0.157594983823183, 12.436944073036383)

    t = TwoSampleBinomialTest(2, 14, 1, 11, 0.1; ptype = :oddratio, htype =:superiority, alpha=0.025,  method=:mover)
    @test_ci_approx t.ci (0.16481005458916503, 15.661396211673583)

    t = TwoSampleBinomialTest(30, 100, 40, 90, 0.1; ptype = :riskratio, htype =:superiority, alpha=0.025,  method=:mn)
    @test_ci_approx t.ci (0.4663609913209758, 0.9799258384796817)

    t = TwoSampleBinomialTest(30, 100, 40, 90, 0.1; ptype = :riskratio, htype =:superiority, alpha=0.025,  method=:cli)
    @test_ci_approx t.ci (0.4663950370893218, 0.9860541079252757)

    t = TwoSampleBinomialTest(30, 100, 40, 90, 0.1; ptype = :riskratio, htype =:superiority, alpha=0.025,  method=:li)
    @test_ci_approx t.ci (0.46246719923578444, 0.9852050064370166)

    t = TwoSampleBinomialTest(30, 100, 40, 90, 0.1; ptype = :riskratio, htype =:superiority, alpha=0.025,  method=:mover)
    @test_ci_approx t.ci (0.46344425873893524, 0.9808806908109405)

    #

    t = TwoSampleBinomialTest(56, 70, 48, 80, 0.1; ptype = :difference, htype = :equivalence, alpha = 0.05,  method = :mn)
    @test_ci_approx t.ci (0.07701994249627694, 0.31666721060091346)

    t = TwoSampleBinomialTest(2, 14, 1, 11, 0.1; ptype = :oddratio, htype = :equivalence, alpha=0.05,  method=:mn)
    @test_ci_approx t.ci (0.2361649043170363, 11.303570631238369)

    t = TwoSampleBinomialTest(30, 100, 40, 90, 0.1; ptype = :riskratio, htype = :equivalence, alpha=0.05,  method=:mn)
    @test_ci_approx t.ci (0.49471139973648776, 0.9228304587889552)

    #

    t = TwoSampleBinomialTest(56, 70, 48, 80, 0.0; ptype = :difference, htype = :equality, alpha = 0.025,  method = :mn)
    @test_ci_approx t.ci (0.031084556392724887, 0.3571127514959887)
    @test t.result == true

    t = TwoSampleBinomialTest(2, 14, 1, 11, 0.0; ptype = :oddratio, htype = :equality, alpha=0.025,  method=:mn)
    @test_ci_approx t.ci (0.13900122510527138, 18.583974295746657)
    @test t.result == true

    t = TwoSampleBinomialTest(30, 100, 40, 90, 0.0; ptype = :riskratio, htype = :equality, alpha=0.025,  method=:mn)
    @test_ci_approx t.ci (0.44251068825683687, 1.0339383423147208)
    @test t.result == true

    #

    show(IOBuffer(), t)
end
