using Base: @deprecate

@deprecate ci(args...) confint(args...)

@deprecate MultinomialLRT MultinomialLRTest
@deprecate OneSampleHotellingT2 OneSampleHotellingT2Test
@deprecate EqualCovHotellingT2 EqualCovHotellingT2Test
@deprecate UnequalCovHotellingT2 UnequalCovHotellingT2Test

@deprecate BartlettsTest BartlettTest

@deprecate confint(x::HypothesisTest, alpha::Real; kwargs...) confint(x; level=1-alpha, kwargs...)