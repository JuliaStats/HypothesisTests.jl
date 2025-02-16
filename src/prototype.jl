using HypothesisTests
using BenchmarkTools

nx = 3
ny = 4

U = 9
@btime HypothesisTests.pwilcox.(0:12, nx, ny, true)
@btime HypothesisTests.pwilcox.((0:12) .-1, nx, ny, false)
##
function _pwilcox(U::Union{Float64, Int}, nx::Int, ny::Int, lower_tail::Bool)
    if !lower_tail
        return _pwilcox(sum((ny+1):(nx+ny)) -sum(1:nx) - U, nx, ny, true)
    end
    Umax = sum((ny+1):(ny+nx))-sum(1:nx)
    if U > Umax÷2
        return 1.0 - _pwilcox(Umax-U-1, nx, ny, true) 
    end
    W = sum((ny + 1):(nx + ny))
    DP = zeros(Int, W, ny + 1)
    for i in 1:(nx+ny)
        for j in i:(ny + 1)
            DP[i,j] = 1
        end
    end
    for k in 2:nx
        for j in 1:(ny+1)
            # for i in U+sum(1:k):-1:1
            for i in min(U+sum(1:k), sum((j):(k+j-1))):-1:1
                if j > 1
                    DP[i, j] = DP[i, j - 1]
                else
                    DP[i, j] = 0
                end
                if i - j - k + 1 >= 1
                    DP[i, j] += DP[i - j - k + 1, j]
                end
            end
        end
    end
    sum(DP[(sum(1:nx)):(sum(1:nx)+U), ny+1]./(prod((ny+1):(nx+ny))÷(prod(1:nx))))
end
@btime _pwilcox.(0:12, nx, ny, true)
@btime _pwilcox.(0:12, nx, ny, false)
