# Dunn test

export DunnTest

struct DunnTest <: HypothesisTest
    df::Int                     # degrees of freedom
    adjustment::Symbol          # adjustment type
    zscores::Vector{Float64}    # z-test statistics
end

"""
    DunnTest(groups::AbstractVector{<:Real}...)

Perform Dunn test for the `groups` ``\\mathcal{G}``. It's a non-parametric pairwise multiple comparisons procedure based on rank sums, often used as a *post hoc* procedure following rejection of a Kruskal–Wallis test. As such, it is a non-parametric analog to multiple pairwise *t* tests following rejection of an ANOVA null hypothesis.

Implements: [`pvalues`](@ref), [`zscores`](@ref)

# References

  * Dunn, O.J., Multiple Comparisons Using Rank Sums. Technometrics, Vol. 6,
    No. 3 (Aug., 1964), pp. 241-252

# External links

  * [Dunn's test on CrossValidated
    ](https://stats.stackexchange.com/tags/dunn-test/info)
"""
function DunnTest(groups::AbstractVector{T}...; adj=:none) where T<:Real
    df = length(groups) - 1
    k = df + 1
    n_i = collect(map(length, groups))
    n = sum(n_i)

    # get ties
    (ranks, tieadj) =  tiedrank_adj([groups...;])
    tieadj /= 12*(length(ranks) - 1)

    # Calculate z-test statistics
    Z = Vector{Float64}(undef, convert(Int,k*(k-1)/2))
    c_i = cumsum(n_i)
    l = 1
    for i in 2:k
        for j in 1:(i-1)
            σ = sqrt( (n*(n+1)/12 - tieadj) * (1/n_i[i] + 1/n_i[j]) )
            mri = mean(view(ranks, (c_i[i-1]+1):c_i[i]))
            mrj = mean(view(ranks, (j == 1 ? 1 : c_i[j-1]+1):c_i[j]))
            Z[l] = (mri - mrj)/σ
            l+=1
        end
    end

    DunnTest(df, adj, Z)
end

testname(::DunnTest) = "Dunn's test"
population_param_of_interest(x::DunnTest) = ("Location parameters", "all equal", NaN)
default_tail(test::DunnTest) = :right
pvalue(t::DunnTest) = NaN
StatsBase.dof(t::DunnTest) = t.df

function show_params(io::IO, x::DunnTest, ident)
    println(io, ident, "Adjustment: ", x.adjustment)
    println(io, ident, "Pairwise comparisons [i - j = Z-score (p-value)]:")
    Z = zscores(x)
    P = pvalues(x)
    for ij in sort!(collect(keys(Z)))
        z = round(Z[ij]; digits=6)
        zpad = z == 0.0 ? " "^5 : ""
        p = round(P[ij]; digits=6)
        ppad = p == 0.0 ? " "^5 : ""
        println(io, ident, "$(ij[1]) - $(ij[2]) = ", (z < 0 ? "$z$zpad" : " $z$zpad"), " ($p$ppad)")
    end
end

# Auxilary methods

"""
    zscores(t::DunnTest)

Returns pairwise z-scores.
"""
zscores(t::DunnTest) = Dict( zip([(j=>i) for i in 2:(dof(t)+1) for j in 1:(i-1)], t.zscores ))

"""
    pvalues(t::DunnTest)

Returns pairwise p-value.
"""
function pvalues(t::DunnTest)
    k = dof(t)+1
    m = k*(k-1)/2
    P = pvalue.(Normal(0.0, 1.0), abs.(t.zscores); tail=default_tail(t))
    Padj = if t.adjustment == :bonferroni
        min.(1, P.*m)
    elseif t.adjustment == :sidak
        min.(1, 1 .- (1 .- P).^m)
    else
        P
    end
    Dict( zip([(j=>i) for i in 2:(dof(t)+1) for j in 1:(i-1)], Padj ))
end
