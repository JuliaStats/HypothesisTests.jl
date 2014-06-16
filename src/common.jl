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

const libRmath = VERSION >= v"0.3.0-" ? "libRmath-julia" : "libRmath"

# PROBABILITY FUNCTIONS
macro libRmath_deferred_free(base)
    libcall = symbol(string(base, "_free"))
    func = symbol(string(base, "_deferred_free"))
    quote
        let gc_tracking_obj = []
            global $func
            function $libcall(x::Vector{None})
                gc_tracking_obj = []
                ccall(($(string(libcall)),libRmath), Void, ())
            end
            function $func()
                if !isa(gc_tracking_obj, Bool)
                    finalizer(gc_tracking_obj, $libcall)
                    gc_tracking_obj = false
                end
            end
        end
    end
end

# RMATH WRAPPERS
macro rmath_deferred_free(base)
    libcall = symbol(string(base, "_free"))
    func = symbol(string(base, "_deferred_free"))
    quote
        let gc_tracking_obj = []
            global $func
            function $libcall(x::Vector{None})
                gc_tracking_obj = []
                ccall(($(string(libcall)),libRmath), Void, ())
            end
            function $func()
                if !isa(gc_tracking_obj, Bool)
                    finalizer(gc_tracking_obj, $libcall)
                    gc_tracking_obj = false
                end
            end
        end
    end
end

## COMMON FUNCTIONS

# Tied rank from Base, modified to compute the adjustment for ties
function tiedrank_adj(v::AbstractArray)
    n     = length(v)
    place = sortperm(v)
    ord   = Array(Float64, n)
    tieadj = 0.0

    i = 1
    while i <= n
        j = i
        while j + 1 <= n && v[place[i]] == v[place[j + 1]]
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