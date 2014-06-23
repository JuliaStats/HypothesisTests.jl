using HypothesisTests, Base.Test

# Seed random generator to avoid random effects in testing
srand(1)

# https://github.com/wojzaremba/btest/blob/master/example.m
X = [8.5415 37.0043; 14.7532 29.2931; 75.5506 28.3312; 29.4645 15.6710; 55.3216 12.3864; 71.9396 43.6335; 32.9615 31.0025; 15.8288 60.9634; 62.5720 61.7679; 59.7134 43.0258; 60.6645 29.5431; 14.0132 45.7089; 58.9875 57.6634; 20.7221 66.7510; 13.4842 56.4697; 59.7120 31.7011; 40.9415 56.8771; 71.6981 71.4109; 30.6346 61.2081; 28.0781 29.2852; 31.0662 30.1430; 57.7878 43.0160; 13.3660 28.4529; 14.4294 43.1398; 60.5189 31.8989; 46.4873 61.3727; 63.1827 16.7833; 29.5423 60.1772; 22.6857 23.5714; 46.9106 44.9033; 32.2429 60.4309; 14.8121 16.1449; 46.5136 14.8556; 12.4426 12.0594; 14.2558 28.6762; 31.3229 76.4065; 8.7645 27.4488; 58.6321 42.0109; 56.5313 70.8958; 32.1743 16.5879]
Y = [75.4855 76.0260; 60.4289 29.7009; 30.0675 59.8129; 60.3929 60.1946; 44.2255 45.7868; 76.4089 29.4659; 74.7562 44.1024; 15.3539 31.5970; 60.9939 15.8542; 59.8324 15.3530; 58.9941 30.7907; 44.0394 58.3662; 76.6321 58.4678; 14.9583 14.3845; 73.2577 15.2053; 58.7344 29.8507; 30.8284 45.2177; 14.6980 31.8136; 76.3094 43.9553; 31.5024 75.7304]
ssize = 20

# test equal distribution
t = BTest(X[1:ssize,:], X[(ssize+1):end,:])
@test t.n == ssize
@test t.block_size == 4
@test_approx_eq pvalue(t) 0.2413024802703047
@test_approx_eq [ci(t)...] [-0.17911833893743978,0.37908317997532043]
show(IOBuffer(), t)

# test different distributions
t = BTest(X[1:ssize,:], Y)
@test_approx_eq pvalue(t) 0.06880259547570035
@test_approx_eq [ci(t)...] [-0.05873950833935862,0.4258075687778003]

# test equal distribution (with laplacian kernel)
t = BTest(X[1:ssize,:], X[(ssize+1):end,:]; kernel=:laplace)
@test_approx_eq pvalue(t) 0.8427600550859411

@test_throws DimensionMismatch BTest(X, Y)
@test_throws ArgumentError BTest(X[1:ssize,:], X[(ssize+1):end,:]; blocksize=1)
@test_throws ArgumentError BTest(X[1:ssize,:], X[(ssize+1):end,:]; kernel=:blubs)