# b.jl
# B-test
#
# Copyright (C) 2014   Christoph Sawade
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

export BTest

immutable BTest <: HypothesisTest
    n::Int           # number of observations
    kernel::Function # kernel function
    block_size::Int  # block size
    MMD::Float64     # distance measure: maximum mean discrepency
    stderr::Float64  # standard error
end

function BTest{T<:Real, S<:Real}(X::AbstractMatrix{T}, Y::AbstractMatrix{S}; kernel::Symbol=:rbf, blocksize::Int=int(floor(sqrt(size(X,1)))))
    n = size(X,1)
    check_blocksize(blocksize, n)
    check_same_size(X, Y)

    if !in(kernel, [:rbf, :laplace, :linear])
        throw(ArgumentError("kernel=$(kernel) is not implemented yet"))
    end

    if kernel == linear
        warn("kernel function is not characteristic")
    end

    kernel_fcn = eval(kernel)
    MMD, stderr = bstats(X,Y, kernel_fcn, blocksize)
    BTest(n, kernel_fcn, blocksize, MMD, stderr)
end

testname(::BTest) = "Block-test"
population_param_of_interest(x::BTest) = ("Maximum mean discrepency", 0, x.MMD) # parameter of interest: name, value under h0, point estimate

function show_params(io::IO, x::BTest, ident)
    println(io, ident, "number of observations: ", x.n)
    println(io, ident, "kernel:                 ", x.kernel)
    println(io, ident, "standard error:         ", x.stderr)
    println(io, ident, "block size:             ", x.block_size)
end

pvalue(x::BTest) = pvalue(Normal(0, x.stderr), x.MMD; tail=:right) # check only if distance is larger than 0
ci(x::BTest, alpha::Float64=0.05) = (quantile(Normal(x.MMD, x.stderr), alpha/2), quantile(Normal(x.MMD, x.stderr), 1-alpha/2))

## helper

# Get block-wise averaged MMD for B-test and its standard deviation
function bstats{T<:Real,S<:Real}(X::AbstractMatrix{T}, Y::AbstractMatrix{S}, kernel::Function, blocksize::Int)
    n = size(X,1)
    X, Y = X[randperm(n),:], Y[randperm(n),:]

    n2::Int = floor(n / blocksize)
    hh = zeros(n2)
    tmp = zeros(n2)
    for i = 1:blocksize
        offset1 = n2*(i-1)
        for j = (i+1):blocksize
            offset2 = n2*(j-1)
            broadcast!(+, hh, hh, kernel(tmp, X, offset1, X, offset2))
            broadcast!(+, hh, hh, kernel(tmp, Y, offset1, Y, offset2))
            broadcast!(-, hh, hh, kernel(tmp, X, offset1, Y, offset2))
            broadcast!(-, hh, hh, kernel(tmp, Y, offset1, X, offset2))
        end
    end
    MMD = mean(hh)
    stderr = std(hh)/sqrt(n2)
    (MMD, stderr)
end

function check_same_size{T<:Real, S<:Real}(X::AbstractMatrix{T}, Y::AbstractMatrix{S})
    if size(X) != size(Y)
        throw(DimensionMismatch("Dimensions of X and Y mismatch."))
    end
end

function check_blocksize(blocksize::Int, n::Int)
    if !(blocksize >= 2 && blocksize < n / 3)
        throw(ArgumentError("block size has to larger than 1 and smaller than $(n/3)"))
    end
end

# point-wise evaluated kernel functions
function l1norm(out::AbstractVector, X::AbstractMatrix, Xoffset::Int, Y::AbstractMatrix, Yoffset::Int)
    fill!(out, 0.)
    for j = 1:size(X, 2), i = 1:length(out)
        out[i] += abs(X[Xoffset+i, j] - Y[Yoffset+i, j])
    end
    out
end
function l2norm2(out::AbstractVector, X::AbstractMatrix, Xoffset::Int, Y::AbstractMatrix, Yoffset::Int)
    fill!(out, 0.)
    for j = 1:size(X, 2), i = 1:length(out)
        out[i] += abs2(X[Xoffset+i, j] - Y[Yoffset+i, j])
    end
    out
end
function rbf(out::AbstractVector, X::AbstractMatrix, Xoffset::Int, Y::AbstractMatrix, Yoffset::Int; sigma::Float64=1.0)
    l2norm2(out, X, Xoffset, Y, Yoffset)
    for i = 1:length(out)
        out[i] = exp(-out[i]/(2*sigma^2))
    end
    out
end
function laplace(out::AbstractVector, X::AbstractMatrix, Xoffset::Int, Y::AbstractMatrix, Yoffset::Int; sigma::Float64=1.0)
    l1norm(out, X, Xoffset, Y, Yoffset)
    for i = 1:length(out)
        out[i] = exp(-out[i]/(2*sigma^2))
    end
    out
end
function linear(out::AbstractVector, X::AbstractMatrix, Xoffset::Int, Y::AbstractMatrix, Yoffset::Int)
    fill!(out, 0.)
    for j = 1:size(X, 2), i = 1:length(out)
        out[i] += X[Xoffset+i, j] * Y[Yoffset+i, j]
    end
    out
end
