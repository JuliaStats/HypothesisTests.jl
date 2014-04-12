using HypothesisTests, Base.Test

@test_approx_eq pvalue(BinomialTest(26, 78)) 0.004334880883507431
@test_approx_eq pvalue(BinomialTest(26, 78), tail=:left) 0.002167440441753716
@test_approx_eq pvalue(BinomialTest(26, 78), tail=:right) 0.9989844298129187
@test_approx_eq [ci(BinomialTest(26, 78))...] [0.2305852396293038,0.4491666887959782]
@test_approx_eq [ci(BinomialTest(26, 78), tail=:left)...] [0.0,0.4313047758370174]
@test_approx_eq [ci(BinomialTest(26, 78), tail=:right)...] [0.2451709633730693,1.0]
show(IOBuffer(), BinomialTest(26, 78))

@test_approx_eq pvalue(BinomialTest([trues(6), falses(3)])) 0.5078125000000002
@test_approx_eq pvalue(BinomialTest([trues(6), falses(3)]), tail=:left) 0.91015625
@test_approx_eq pvalue(BinomialTest([trues(6), falses(3)]), tail=:right) 0.2539062500000001
@test_approx_eq [ci(BinomialTest([trues(6), falses(3)]))...] [0.2992950562085405,0.9251453685803082]
@test_approx_eq [ci(BinomialTest([trues(6), falses(3)]), tail=:left)...] [0.0,0.9022531865607242]
@test_approx_eq [ci(BinomialTest([trues(6), falses(3)]), tail=:right)...] [0.3449413659437032,1.0]
show(IOBuffer(), BinomialTest([trues(6), falses(3)]))

x = [55, 58, 61, 61, 62, 62, 62, 63, 63, 64, 66, 68, 68, 69, 69, 69, 70, 71, 72, 72]
@test_approx_eq pvalue(SignTest(x)) 1.907348632812499e-6
@test_approx_eq pvalue(SignTest(x), tail=:left) 1.0
@test_approx_eq pvalue(SignTest(x), tail=:right) 9.536743164062495e-7
@test_approx_eq pvalue(SignTest(x, 70)) 0.004425048828125003
@test_approx_eq pvalue(SignTest(x, 70), tail=:left) 0.0022125244140625013
@test_approx_eq pvalue(SignTest(x, 70), tail=:right) 0.9996356964111328
show(IOBuffer(), SignTest(x, 70))

x = [9, 2, 7, 5]
y = [7, 2, 6, 4]
@test_approx_eq pvalue(SignTest(x, y)) 0.25
show(IOBuffer(), SignTest(x, y))
