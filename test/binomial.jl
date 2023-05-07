using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "Binomial" begin
    t = BinomialTest(26, 78)
    @test pvalue(t) ≈ 0.004334880883507431
    @test pvalue(t, tail=:left) ≈ 0.002167440441753716
    @test pvalue(t, tail=:right) ≈ 0.9989844298129187
    @test default_tail(t) == :both
    @test_ci_approx confint(t) (0.23058523962930383, 0.4491666887959782)
    @test_ci_approx confint(t, tail=:left) (0.0, 0.4313047758370174)
    @test_ci_approx confint(t, tail=:right) (0.2451709633730693, 1.0)
    @test_ci_approx confint(t, method=:wald) (0.22871819521037956, 0.43794847145628707)
    @test_ci_approx confint(t, tail=:left, method=:wald) (0.0, 0.42112912485444692)
    @test_ci_approx confint(t, tail=:right, method=:wald) (0.24553754181221971, 1.0)
    @test_ci_approx confint(t, method=:wilson) (0.23872670036358601, 0.44358590287381217)
    @test_ci_approx confint(t, tail=:left, method=:wilson) (0.0, 0.42541288951088108)
    @test_ci_approx confint(t, tail=:right, method=:wilson) (0.25242832328277831, 1.0)
    @test_ci_approx confint(t, method=:jeffrey) (0.23626570247518358, 0.44251318323879296)
    @test_ci_approx confint(t, tail=:left, method=:jeffrey) (0.0, 0.42466492683653623)
    @test_ci_approx confint(t, tail=:right, method=:jeffrey) (0.25098836986261724, 1.0)
    @test_ci_approx confint(t, method=:agresti_coull) (0.2384423809121706, 0.44387022232522744)
    @test_ci_approx confint(t, tail=:left, method=:agresti_coull) (0.0, 0.42558712894362222)
    @test_ci_approx confint(t, tail=:right, method=:agresti_coull) (0.25225408385003706, 1.0)
    @test_ci_approx confint(t, method=:arcsine) (0.23366209634204066,0.44117918327686334)
    @test_ci_approx confint(t, tail=:left, method=:arcsine) (0.0,0.4235046425920888)
    @test_ci_approx confint(t, tail=:right, method=:arcsine) (0.2489264087216164,1.0)

    show(IOBuffer(), t)

    t = BinomialTest([trues(6); falses(3)])
    @test pvalue(t) ≈ 0.5078125000000002
    @test pvalue(t, tail=:left) ≈ 0.91015625
    @test pvalue(t, tail=:right) ≈ 0.2539062500000001
    @test_ci_approx confint(t) (0.2992950562085405, 0.9251453685803082)
    @test_ci_approx confint(t, tail=:left) (0.0, 0.9022531865607242)
    @test_ci_approx confint(t, tail=:right) (0.3449413659437032, 1.0)
    @test_ci_approx confint(t, method=:wald) (0.35868803903340479, 0.97464529429992841)
    @test_ci_approx confint(t, tail=:left, method=:wald) (0.0,0.92513047859481645)
    @test_ci_approx confint(t, tail=:right, method=:wald) (0.40820285473851681,1.0)
    @test_ci_approx confint(t, method=:wilson) (0.35420213558039609,0.87941618161308899)
    @test_ci_approx confint(t, tail=:left, method=:wilson) (0.0,0.85802909820500495)
    @test_ci_approx confint(t, tail=:right, method=:wilson) (0.39825972868840931,1.0)
    @test_ci_approx confint(t, method=:jeffrey) (0.34779179347226591,0.89578677833922582)
    @test_ci_approx confint(t, tail=:left, method=:jeffrey) (0.0,0.86830830610561005)
    @test_ci_approx confint(t, tail=:right, method=:jeffrey) (0.39604343455469687,1.0)
    @test_ci_approx confint(t, method=:agresti_coull) (0.350905767251112,0.88271254994237336)
    @test_ci_approx confint(t, tail=:left, method=:agresti_coull) (0.0,0.86049746046629294)
    @test_ci_approx confint(t, tail=:right, method=:agresti_coull) (0.39579136642712159,1.0)
    @test_ci_approx confint(t, method=:arcsine) (0.345812446615087,0.9188773496172281)
    @test_ci_approx confint(t, tail=:left, method=:arcsine) (0.0,0.8879439981269358)
    @test_ci_approx confint(t, tail=:right, method=:arcsine) (0.3965293068864491,1.0)
    show(IOBuffer(), t)

    t = BinomialTest(0, 100, 0.01)
    @test pvalue(t) ≈ 0.7320646825464591
    show(IOBuffer(), t)

    t = BinomialTest(100, 100, 0.99)
    @test pvalue(t) ≈ 0.7320646825464584
    show(IOBuffer(), t)

    # from issue #295
    # without clamping: (-0.05457239484968546, 0.4890548596328611)
    @test_ci_approx confint(BinomialTest(0, 5), method=:agresti_coull) (0.0, 0.4890548596328611)
    # without clamping: (0.5109451403671388, 1.0545723948496855)
    @test_ci_approx confint(BinomialTest(5, 5), method=:agresti_coull) (0.5109451403671388, 1.0)
    # without clamping: (-0.15060901623063327, 0.5506090162306333)
    @test_ci_approx confint(BinomialTest(1, 5), method=:wald) (0.0, 0.5506090162306333)
    # without clamping: (0.44939098376936687, 1.1506090162306333)
    @test_ci_approx confint(BinomialTest(4, 5), method=:wald) (0.44939098376936687, 1.0)
    # without clamping: (-2.7755575615628914e-17, 0.2775327998628899)
    @test_ci_approx confint(BinomialTest(0, 10), method=:wilson) (0.0, 0.2775327998628899)
    # without clamping: (0.7575059933447587, 1.0000000000000002)
    @test_ci_approx confint(BinomialTest(12, 12), method=:wilson) (0.7575059933447587, 1.0)
end

@testset "SignTest" begin
    x = [55, 58, 61, 61, 62, 62, 62, 63, 63, 64, 66, 68, 68, 69, 69, 69, 70, 71, 72, 72]
    @test pvalue(SignTest(x)) ≈ 1.907348632812499e-6
    @test pvalue(SignTest(x), tail=:left) ≈ 1.0
    @test pvalue(SignTest(x), tail=:right) ≈ 9.536743164062495e-7
    @test pvalue(SignTest(x, 70)) ≈ 0.004425048828125003
    @test pvalue(SignTest(x, 70), tail=:left) ≈ 0.0022125244140625013
    @test pvalue(SignTest(x, 70), tail=:right) ≈ 0.9996356964111328
    @test default_tail(SignTest(x)) == :both
    @test_ci_approx confint(SignTest(x, 70)) (62, 69)
    @test_ci_approx confint(SignTest(x, 70), level=0.9998) (61, 71)
    show(IOBuffer(), SignTest(x, 70))

    x = [9, 2, 7, 5]
    y = [7, 2, 6, 4]
    @test pvalue(SignTest(x, y)) ≈ 0.25
    @test_ci_approx confint(SignTest(x, y)) (0, 2)
    show(IOBuffer(), SignTest(x, y))

    # www.stat.umn.edu/geyer/old03/5102/notes/rank.pdf
    x = [-4.7, 3.7, 22.4, 13.6, 8.7, 9.1, -7.8, 10.8, 15.6, 23.5, 14.4, 20.2, 6.5, 10.1, -6.9]
    @test pvalue(SignTest(x)) ≈ 0.03515625
    @test pvalue(SignTest(x), tail=:left) ≈ 0.996307373046875
    @test pvalue(SignTest(x), tail=:right) ≈ 0.017578125000000007
    @test_ci_approx confint(SignTest(x), level=0.99995) (-7.8, 23.5)
    @test_ci_approx confint(SignTest(x), level=0.9999) (-6.9, 22.4)
    @test_ci_approx confint(SignTest(x), level=0.993) (-4.7, 20.2)
    @test_ci_approx confint(SignTest(x), level=0.965) (3.7, 15.6)
    @test_ci_approx confint(SignTest(x), level=0.882) (6.5, 14.4)
    @test_ci_approx confint(SignTest(x), level=0.7) (8.7, 13.6)
    @test_ci_approx confint(SignTest(x), level=0.6) (9.1, 10.8)
    show(IOBuffer(), SignTest(x))
end
