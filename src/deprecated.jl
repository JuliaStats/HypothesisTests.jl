using Base: @deprecate

@deprecate ci(args...) confint(args...)

@deprecate MultinomialLRT MultinomialLRTest
@deprecate OneSampleHotellingT2 OneSampleHotellingT2Test
@deprecate EqualCovHotellingT2 EqualCovHotellingT2Test
@deprecate UnequalCovHotellingT2 UnequalCovHotellingT2Test