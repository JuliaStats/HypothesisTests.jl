load("test")
load("src/HypothesisTests")
using HypothesisTests, Test

# Fisher example 4.11
@test abs(p_value(RayleighTest, 60*0.2370^2, 60) - 0.034) <= 0.001
# Fisher example 4.12
@test abs(p_value(RayleighTest, [2, 9, 18, 24, 30, 35, 35, 39, 39, 44,
 44, 49, 56, 70, 76, 76, 81, 86, 91, 112,
 121, 127, 133, 134, 138, 147, 152, 157, 166, 171,
 177, 187, 206, 210, 211, 215, 238, 246, 269, 270,
 285, 292, 305, 315, 325, 328, 329, 343, 354, 359]*pi/180) - 0.20) <= 0.01