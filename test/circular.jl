load("test")
load("src/HypothesisTests")
using HypothesisTests, Test

# Fisher, 1995 example 4.11
@test abs(p_value(RayleighTest, 60*0.2370^2, 60) - 0.034) <= 0.001

# Fisher, 1995 example 4.12
@test abs(p_value(RayleighTest,
    [2, 9, 18, 24, 30, 35, 35, 39, 39, 44,
	44, 49, 56, 70, 76, 76, 81, 86, 91, 112,
	121, 127, 133, 134, 138, 147, 152, 157, 166, 171,
	177, 187, 206, 210, 211, 215, 238, 246, 269, 270,
	285, 292, 305, 315, 325, 328, 329, 343, 354, 359]
	*pi/180) - 0.20) <= 0.01

# Fisher, 1995 example 6.8
wind_direction_6am =
	[356, 97, 211, 232, 343, 292, 157, 302, 335, 302,
	324, 85, 324, 340, 157, 238, 254, 146, 232, 122,
	329]*pi/180
wind_direction_12pm =
	[119, 162, 221, 259, 270, 29, 97, 292, 40, 313,
	94, 45, 47, 108, 221, 270, 119, 248, 270, 45,
	23]*pi/180
@test abs(test_statistic(FisherTLinearAssociation, wind_direction_6am, wind_direction_12pm) - 0.191) < 0.001
@test abs(p_value(FisherTLinearAssociation, wind_direction_6am, wind_direction_12pm) - 0.01) < 0.01

# Jammaladak, 2001 example 8.1
@test abs(test_statistic(JammaladakCircularCorrelation, wind_direction_6am, wind_direction_12pm) - 0.2704648) < 1e-7
@test abs(p_value(JammaladakCircularCorrelation, wind_direction_6am, wind_direction_12pm) - 0.2247383) < 1e-7