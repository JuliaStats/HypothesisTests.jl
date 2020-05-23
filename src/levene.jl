# Levene Test for Equality of Variances

export LeveneTest

abstract type LeveneTest <: HypothesisTest end

population_param_of_interest(::LeveneTest) =
    ("Equality of variances", NaN, NaN)

