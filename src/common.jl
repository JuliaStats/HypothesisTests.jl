# common.jl
# Common provides types and helper functions
#
# Copyright (C) 2012   Simon Kornblith
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


## COMMON FUNCTIONS

# Tied ranking, also computing the adjustment for ties
# sum(n^3-n) where n is the number of ties for each data point
function tiedrank_adj!(ord::AbstractVector, v::AbstractArray)
    n     = length(v)
    (n == length(ord)) ||
        throw(DimensionMismatch("The length of ord ($(length(ord))) differs from v length ($n)"))
    place = sortperm(v)
    tieadj = 0.0

    i = 1
    @inbounds while i <= n
        j = i
        vi = v[place[i]]
        while (j + 1 <= n) && (vi == v[place[j + 1]])
            j += 1
        end

        if j > i
            t = j - i + 1
            m = sum(i:j) / t
            tieadj += t^3 - t
            for k = i:j
                ord[place[k]] = m
            end
        else
            ord[place[i]] = i
        end

        i = j + 1
    end

    (ord, tieadj)
end

tiedrank_adj(v::AbstractArray) = tiedrank_adj!(Vector{Float64}(undef, length(v)), v)

# Pool covariance matrices, overwriting the first
function poolcov!(Sx::AbstractMatrix, nxm1::Int, Sy::AbstractMatrix, nym1::Int)
    @inbounds for i = eachindex(Sx, Sy)
        Sx[i] = (Sx[i] * nxm1 + Sy[i] * nym1) / (nxm1 + nym1)
    end
    Sx
end
