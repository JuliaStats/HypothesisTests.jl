using Base: @deprecate

@deprecate ci(args...) confint(args...)
@deprecate default_tail(test::HypothesisTest) tail(test)
