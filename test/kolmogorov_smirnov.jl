using HypothesisTests, Distributions, Base.Test

# sample drawn from uniform distribution
x = [0.3500, 0.1966, 0.2511, 0.6160, 0.4733,
     0.3517, 0.8308, 0.5853, 0.5497, 0.9172,
     0.2858, 0.7572, 0.7537, 0.3804, 0.5678,
     0.0759, 0.0540, 0.5308, 0.7792, 0.9340,
     0.1299, 0.5688, 0.4694, 0.0119, 0.3371
]
t = ApproximateOneSampleKSTest(x, Uniform())
@test_approx_eq t.δ 0.1440
@test_approx_eq t.δn 0.0571
@test_approx_eq t.δp 0.1440
@test_approx_eq pvalue(t) 0.6777349664784745
@test_approx_eq pvalue(t; tail=:left) 0.849573771973747
@test_approx_eq pvalue(t; tail=:right) 0.3545875485608989
show(IOBuffer(), t)

t = ApproximateTwoSampleKSTest(x, [(0:24)/25...])
@test_approx_eq t.δ 0.12
@test_approx_eq t.δn 0.08
@test_approx_eq t.δp 0.12
@test_approx_eq pvalue(t) 0.993764859699076
@test_approx_eq pvalue(t; tail=:left) 0.8521437889662113
@test_approx_eq pvalue(t; tail=:right) 0.697676326071031
show(IOBuffer(), t)

t = ExactOneSampleKSTest(x, Uniform())
@test_approx_eq t.δ 0.1440
@test_approx_eq t.δn 0.0571
@test_approx_eq t.δp 0.1440
@test_approx_eq pvalue(t) 0.6263437768244742
@test_approx_eq pvalue(t; tail=:left) 0.8195705417998183
@test_approx_eq pvalue(t; tail=:right) 0.32350648882777194
show(IOBuffer(), t)

## check fit to normal distribution
t = ApproximateOneSampleKSTest(x, Normal())
@test_approx_eq t.δ 0.5047473010922947
@test_approx_eq t.δn 0.5047473010922947
@test_approx_eq t.δp 0.17515194649718513
@test_approx_eq pvalue(t) 5.871827067532435e-6
@test_approx_eq pvalue(t; tail=:left) 2.9359135337662175e-6
@test_approx_eq pvalue(t; tail=:right) 0.21569061887162347

## check unequal sample size
t = ApproximateTwoSampleKSTest(x, [(0:5)/6...])
@test_approx_eq t.δ 0.22
@test_approx_eq t.δn 0.22
@test_approx_eq t.δp 0.09333333333333346
@test_approx_eq pvalue(t) 0.973300892518972
@test_approx_eq pvalue(t; tail=:left) 0.6260111498528065
@test_approx_eq pvalue(t; tail=:right) 0.9191544797498837

# http://ocw.mit.edu/courses/mathematics/18-443-statistics-for-applications-fall-2006/lecture-notes/lecture14.pdf 
x = [0.58, 0.42, 0.52, 0.33, 0.43, 0.23, 0.58, 0.76, 0.53, 0.64]
t = ApproximateOneSampleKSTest(x, Uniform())
@test_approx_eq t.δ 0.26
@test_approx_eq t.δn 0.23
@test_approx_eq t.δp 0.26
@test_approx_eq pvalue(t) 0.5084937988981307
@test_approx_eq pvalue(t; tail=:left) 0.3471494153245104
@test_approx_eq pvalue(t; tail=:right) 0.25872229825964005

t = ApproximateTwoSampleKSTest(x, [(0:9)/10...])
@test_approx_eq t.δ 0.3
@test_approx_eq t.δn 0.3
@test_approx_eq t.δp 0.2
@test_approx_eq pvalue(t) 0.7590978384203948
@test_approx_eq pvalue(t; tail=:left) 0.406569659740599
@test_approx_eq pvalue(t; tail=:right) 0.6703200460356393

t = ExactOneSampleKSTest(x, Uniform())
@test_approx_eq t.δ 0.26
@test_approx_eq t.δn 0.23
@test_approx_eq t.δp 0.26
@test_approx_eq pvalue(t) 0.4351284228580825
@test_approx_eq pvalue(t; tail=:left) 0.3013310572470338
@test_approx_eq pvalue(t; tail=:right) 0.2193143479950862
