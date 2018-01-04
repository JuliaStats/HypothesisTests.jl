#=
From:
PATRICK ROYSTON
Approximating the Shapiro-Wilk W-test for non-normality
*Statistics and Computing* (1992) **2**, 117-119
DOI: [10.1007/BF01891203](https://doi.org/10.1007/BF01891203)
=#

# TODO: Rerun simulation and polynomial fitting

const ROYSTON_COEFFS = Dict{String, Vector{Float64}}(
"C1" => [0.0E0, 0.221157E0, -0.147981E0, -0.207119E1, 0.4434685E1, -0.2706056E1],
"C2" => [0.0E0, 0.42981E-1, -0.293762E0, -0.1752461E1, 0.5682633E1, -0.3582633E1],
"C3" => [0.5440E0, -0.39978E0, 0.25054E-1, -0.6714E-3],
"C4" => [0.13822E1, -0.77857E0, 0.62767E-1, -0.20322E-2],
"C5" => [-0.15861E1, -0.31082E0, -0.83751E-1, 0.38915E-2],
"C6" => [-0.4803E0, -0.82676E-1, 0.30302E-2],
"C7" => [0.164E0, 0.533E0],
"C8" => [0.1736E0, 0.315E0],
"C9" => [0.256E0, -0.635E-2],
"G"  => [-0.2273E1, 0.459E0]
)

for (s,c) in ROYSTON_COEFFS
    @eval $(Symbol("_"*s))(x) = Base.Math.@horner(x, $(c...))
end

#=
The following hardcoded constants has been replaced by more precise values:

SQRTH = sqrt(2.0)/2.0 # 0.70711E0
TH = 3/8 # 0.375E0
SMALL = eps(1.0) # 1E-19
PI6 = π/6 # 0.1909859E1
STQR = asin(sqrt(0.75)) # 0.1047198E1
#=
