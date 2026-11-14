using Base.Test

my_tests = (
    "binomial.jl",
    "kolmogorov_smirnov.jl",
    "kruskal_wallis.jl",
    "mann_whitney.jl",
    "t.jl",
    "z.jl",
    # "wilcoxon.jl",
    "anderson_darling.jl",
)

println("Running simulation tests:")

@testset "All simulation studies" begin
    for my_test in my_tests
        try
            include(my_test)
            println("\t\033[1m\033[32mPASSED\033[0m: $(my_test)")
        catch e
            println("\t\033[1m\033[31mFAILED\033[0m: $(my_test)")
        end
    end
end
