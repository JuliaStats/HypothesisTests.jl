# The t-tests

Specification of three procedures for inference on means:

  - the **one-sample t-test**, for a single sample against a hypothesised mean, and its
    **paired** form, which is the one-sample test applied to within-pair differences;
  - the **two-sample t-test assuming equal variances**, Student's;
  - the **two-sample t-test not assuming equal variances**, Welch's.

All three share their p-value ([§6](@ref "6. p-values")) and their confidence interval
([§7](@ref "7. Confidence interval")); they differ only in what the estimate, the standard
error, and the degrees of freedom are
([§3](@ref "3. One-sample")–[§5](@ref "5. Two samples, unequal variances (Welch)")).
[§8](@ref "8. Worked values for the t-tests") gives worked values for checking an implementation.

**In this package.** The one-sample test, and its paired form, is
[`OneSampleTTest`](@ref); the two-sample tests are [`EqualVarianceTTest`](@ref) and
[`UnequalVarianceTTest`](@ref). [§9](@ref "9. The t-tests in this package") maps the specification
onto them and records where they depart from it.

## 1. Preliminaries and notation

| symbol | meaning |
|:---|:---|
| ``x_1, \dots, x_{n}`` | the one-sample input; for paired data these are the differences ``x_i - y_i`` |
| ``x_1, \dots, x_{n_x}``, ``y_1, \dots, y_{n_y}`` | the two-sample inputs |
| ``\bar x`` | the sample mean ``\frac{1}{n}\sum_i x_i`` |
| ``s^2`` | the unbiased sample variance ``\frac{1}{n-1}\sum_i (x_i - \bar x)^2`` |
| ``\mu_0`` | the value of the estimand under the null hypothesis |
| ``\hat\delta`` | the point estimate: ``\bar x``, or ``\bar x - \bar y`` |
| ``\mathrm{SE}`` | the estimated standard error of ``\hat\delta`` |
| ``\nu`` | degrees of freedom |
| ``\alpha`` | two-sided error rate; the interval has coverage ``1-\alpha`` |
| ``T_\nu`` | Student's ``t`` distribution on ``\nu`` degrees of freedom |
| ``t_{\nu, q}`` | its ``q`` quantile |

Note that ``s^2`` carries the ``n-1`` denominator throughout. The sample variance is
sometimes defined with an ``n`` denominator instead; an implementation that uses that
form will disagree with every standard error, statistic, p-value and interval on this
page, which makes it the first thing to check when conformance fails.

## 2. Model and estimand

**Model.** The observations are independent and normally distributed with a common mean
and a common variance within each sample. Only the mean is being tested; the variance is a
nuisance parameter, estimated from the data, and it is that estimation which produces a
``t`` rather than a normal reference distribution.

**Estimand.** The population mean ``\mu`` for the one-sample test, and the difference of
population means ``\mu_x - \mu_y`` for the two-sample tests.

**Null hypothesis.** That the estimand equals ``\mu_0``, which defaults to ``0``.

**What each test assumes.**

| test | assumes |
|:---|:---|
| one-sample ([§3](@ref "3. One-sample")) | normality of the ``x_i`` |
| paired ([§3.1](@ref "3.1 Paired")) | normality of the *differences*, not of either sample |
| Student ([§4](@ref "4. Two samples, equal variances (Student)")) | normality of both samples, and ``\sigma_x^2 = \sigma_y^2`` |
| Welch ([§5](@ref "5. Two samples, unequal variances (Welch)")) | normality of both samples |

Normality matters less than it appears: by the central limit theorem the distribution of
``\hat\delta`` approaches normality as the sample grows, and the ``t`` reference
distribution approaches the normal, so the test is asymptotically valid for any
distribution with finite variance. What does not wash out is dependence between
observations, and, for [§4](@ref "4. Two samples, equal variances (Student)"), a
difference in variances when the samples are also of different sizes.

## 3. One-sample

```math
\hat\delta = \bar x , \qquad
\mathrm{SE} = \frac{s}{\sqrt{n}} , \qquad
\nu = n - 1 , \qquad
t = \frac{\hat\delta - \mu_0}{\mathrm{SE}} .
```

Under the null, ``t \sim T_{n-1}`` exactly when the ``x_i`` are normal.

### 3.1 Paired

The paired test is not a separate procedure: it is [§3](@ref "3. One-sample") applied to
the differences ``d_i = x_i - y_i``, which requires the two inputs to be of equal length
and in corresponding order. Since ``\overline{x-y} = \bar x - \bar y``, the estimate
agrees with the two-sample one; the standard error does not, because pairing removes the
between-subject variance from it. Pairing is a property of the data.

## 4. Two samples, equal variances (Student)

The two samples are assumed to share a variance, so both contribute to a single pooled
estimate of it, weighted by degrees of freedom:

```math
s_p^2 = \frac{(n_x - 1)s_x^2 + (n_y - 1)s_y^2}{n_x + n_y - 2} .
```

```math
\hat\delta = \bar x - \bar y , \qquad
\mathrm{SE} = s_p \sqrt{\frac{1}{n_x} + \frac{1}{n_y}} , \qquad
\nu = n_x + n_y - 2 , \qquad
t = \frac{\hat\delta - \mu_0}{\mathrm{SE}} .
```

Under the null and the stated assumptions, ``t \sim T_{n_x + n_y - 2}`` exactly.

## 5. Two samples, unequal variances (Welch)

Each variance is estimated from its own sample, and no pooling occurs:

```math
\hat\delta = \bar x - \bar y , \qquad
\mathrm{SE} = \sqrt{\frac{s_x^2}{n_x} + \frac{s_y^2}{n_y}} , \qquad
t = \frac{\hat\delta - \mu_0}{\mathrm{SE}} .
```

The price is that ``t`` is no longer exactly ``t``-distributed under any ``\nu``. The
Welch-Satterthwaite equation supplies the ``\nu`` that matches the first two moments of the
variance estimator:

```math
\nu = \frac{\left(\dfrac{s_x^2}{n_x} + \dfrac{s_y^2}{n_y}\right)^{\!2}}
           {\dfrac{(s_x^2/n_x)^2}{n_x - 1} + \dfrac{(s_y^2/n_y)^2}{n_y - 1}} .
```

``\nu`` is not in general an integer and must not be rounded; ``T_\nu`` is defined for real
``\nu > 0``. It satisfies ``\min(n_x, n_y) - 1 \le \nu \le n_x + n_y - 2``, reaching the
upper bound when the two variance estimates and sample sizes coincide.

## 6. p-values

Identical for all three, with ``\nu`` and ``t`` from the relevant section. Writing ``F``
for the CDF of ``T_\nu``,

```math
p_{\text{left}} = F(t), \qquad
p_{\text{right}} = 1 - F(t), \qquad
p_{\text{both}} = 2\bigl(1 - F(|t|)\bigr) .
```

The two-sided form uses the symmetry of ``T_\nu`` about zero and needs no clipping, since
``F(|t|) \ge 1/2``. `left` is the alternative that the estimand lies below ``\mu_0``, and
`right` that it lies above.

## 7. Confidence interval

By inversion of the two-sided test: the set of ``\mu_0`` not rejected at level ``\alpha``.
Since ``t`` is monotone in ``\mu_0``, that set is the interval

```math
\hat\delta \pm t_{\nu,\, 1 - \alpha/2} \cdot \mathrm{SE} ,
```

symmetric about the estimate and independent of ``\mu_0``. Its coverage is exact under the
model.

### 7.1 One-sided intervals

A one-sided bound at level ``L`` is the corresponding endpoint of the two-sided interval at
level ``2L - 1``, the other endpoint being infinite. Which endpoint is retained follows the
alternative the tail names: the alternative ``\mu < \mu_0`` is compatible with an upper
bound, so

```math
\text{left} \;\longrightarrow\; \bigl(-\infty,\; \hat\delta + t_{\nu, L}\,\mathrm{SE}\bigr) ,
\qquad
\text{right} \;\longrightarrow\; \bigl(\hat\delta - t_{\nu, L}\,\mathrm{SE},\; \infty\bigr) .
```

This is the convention R uses, under `alternative = "less"` and `"greater"`, and the one
every test in this package follows, the rank tests included: see
[§6.5](@ref "6.5 One-sided intervals") of [Rank-based location inference](@ref).

## 8. Worked values for the t-tests

Conformance vectors; the table values are printed to six decimal places. The sessions are run when this page is
built, so what is shown is what the package returns; the tables give the intermediate
quantities that no printed output shows. One pair of samples serves the first three tests:

``x`` = `[5.1, 4.9, 6.2, 5.8, 5.3, 6.1, 5.5, 5.9, 4.7, 6.0]`, ``n_x = 10``;
``y`` = `[4.8, 5.2, 4.5, 5.0, 4.9, 5.4, 4.6, 5.1]`, ``n_y = 8``.

The paired test needs a partner of equal length, so it gets its own pair below.

```jldoctest ttest
julia> using HypothesisTests

julia> x = [5.1, 4.9, 6.2, 5.8, 5.3, 6.1, 5.5, 5.9, 4.7, 6.0];

julia> y = [4.8, 5.2, 4.5, 5.0, 4.9, 5.4, 4.6, 5.1];
```

**One sample against ``\mu_0 = 5``** ([§3](@ref "3. One-sample")). ``\bar x = 5.55``,
``s = 0.529675``. The null value is the trailing positional argument, so testing a
different ``\mu_0`` is a new call on the same data.

```jldoctest ttest
julia> t = OneSampleTTest(x, 5)
One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         5
    point estimate:          5.55
    95% confidence interval: (5.171, 5.929)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0095

Details:
    number of observations:   10
    t-statistic:              3.2836227276929635
    degrees of freedom:       9
    empirical standard error: 0.1674979270186815


julia> pvalue(t)
0.00947430516135569

julia> pvalue(t; tail = :right)
0.004737152580677845

julia> confint(t)
(5.171093364640838, 5.928906635359161)

julia> confint(t; level = 0.90)
(5.242957383788944, 5.857042616211055)
```

| quantity | value |
|:---|:---|
| ``\mathrm{SE}`` | `0.167498` |
| ``\nu`` | `9` |
| ``t`` | `3.283623` |
| ``p_{\text{both}}`` | `0.009474` |
| ``p_{\text{right}}`` | `0.004737` |
| interval, ``1-\alpha = 0.95`` | `(5.171093, 5.928907)` |
| interval, ``1-\alpha = 0.90`` | `(5.242957, 5.857043)` |

with ``t_{9,\,0.975} = 2.262157``. Note ``p_{\text{both}} = 2 p_{\text{right}}`` exactly,
which [§6](@ref "6. p-values") gets from the symmetry of the ``t`` distribution and no
clipping: unlike the rank tests, the statistic here is continuous, so a doubled tail cannot
exceed ``1``.

**Two samples, equal variances** ([§4](@ref "4. Two samples, equal variances (Student)")).
``\hat\delta = 0.6125``, ``s_p = 0.444673``.

```jldoctest ttest
julia> t = EqualVarianceTTest(x, y)
Two sample t-test (equal variance)
----------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          0.6125
    95% confidence interval: (0.1654, 1.06)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0104

Details:
    number of observations:   [10,8]
    t-statistic:              2.9038471071017558
    degrees of freedom:       16
    empirical standard error: 0.210927083076119


julia> confint(t)
(0.1653545588376535, 1.0596454411623462)
```

| quantity | value |
|:---|:---|
| ``\mathrm{SE}`` | `0.210927` |
| ``\nu`` | `16` |
| ``t`` | `2.903847` |
| ``p_{\text{both}}`` | `0.010358` |
| interval, ``1-\alpha = 0.95`` | `(0.165355, 1.059645)` |

**Two samples, unequal variances**
([§5](@ref "5. Two samples, unequal variances (Welch)")). ``\hat\delta = 0.6125``. Same
data, and it must be asked for by name, since nothing dispatches to it.

```jldoctest ttest
julia> t = UnequalVarianceTTest(x, y)
Two sample t-test (unequal variance)
------------------------------------
Population details:
    parameter of interest:   Mean difference
    value under h_0:         0
    point estimate:          0.6125
    95% confidence interval: (0.1883, 1.037)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0077

Details:
    number of observations:   [10,8]
    t-statistic:              3.0833130203886516
    degrees of freedom:       14.684901529359335
    empirical standard error: 0.19864995735100363


julia> confint(t)
(0.1882949619174023, 1.0367050380825973)
```

| quantity | value |
|:---|:---|
| ``\mathrm{SE}`` | `0.198650` |
| ``\nu`` | `14.684902` |
| ``t`` | `3.083313` |
| ``p_{\text{both}}`` | `0.007727` |
| interval, ``1-\alpha = 0.95`` | `(0.188295, 1.036705)` |

Note ``\nu`` between ``\min(n_x,n_y) - 1 = 7`` and ``n_x + n_y - 2 = 16``, and the Welch
interval narrower here than the Student one, because the larger sample carries the larger
variance.

**Paired** ([§3.1](@ref "3.1 Paired")). The same ``x``, now paired with a second
ten-point sample ``y`` = `[5.0, 4.6, 6.0, 5.5, 5.1, 5.8, 5.2, 5.6, 4.4, 5.7]`, so the
differences of [§3.1](@ref "3.1 Paired") are ``d_i = x_i - y_i``. The two-argument
`OneSampleTTest` is the paired form: it computes ``d`` itself, and calling it on the
differences is the same test. This is a session of its own, since ``y`` here is a
different sample from the ``y`` above.

```jldoctest ttestpaired
julia> using HypothesisTests

julia> x = [5.1, 4.9, 6.2, 5.8, 5.3, 6.1, 5.5, 5.9, 4.7, 6.0];

julia> y = [5.0, 4.6, 6.0, 5.5, 5.1, 5.8, 5.2, 5.6, 4.4, 5.7];

julia> d = x .- y;

julia> t = OneSampleTTest(x, y)
One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0
    point estimate:          0.26
    95% confidence interval: (0.21, 0.31)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           <1e-06

Details:
    number of observations:   10
    t-statistic:              11.75894243853277
    degrees of freedom:       9
    empirical standard error: 0.022110831935702693


julia> confint(t)
(0.20998182316122294, 0.3100181768387772)

julia> pvalue(OneSampleTTest(d)) == pvalue(t)
true
```

| quantity | value |
|:---|:---|
| ``\hat\delta`` | `0.260000` |
| ``\nu`` | `9` |
| ``t`` | `11.758942` |
| interval, ``1-\alpha = 0.95`` | `(0.209982, 0.310018)` |

The standard error here, `0.022111`, is far below the one-sample figure above, because
pairing removes the between-subject variation: that is the point of
[§3.1](@ref "3.1 Paired").

## 9. The t-tests in this package

[§3](@ref "3. One-sample") is [`OneSampleTTest`](@ref), whose paired form
[§3.1](@ref "3.1 Paired") is the two-argument method.
[§4](@ref "4. Two samples, equal variances (Student)") is [`EqualVarianceTTest`](@ref) and
[§5](@ref "5. Two samples, unequal variances (Welch)") is [`UnequalVarianceTTest`](@ref);
both are subtypes of `TwoSampleTTest`. Unlike R, whose `t.test` is a single entry point
with a `var.equal` switch defaulting to Welch's, this package has no entry point that
chooses between them, so [§5](@ref "5. Two samples, unequal variances (Welch)")
must be asked for by name.

`pvalue` implements [§6](@ref "6. p-values") and `confint` implements
[§7](@ref "7. Confidence interval"), both on the shared `TTest` supertype. All three tests
accept `μ0` as a trailing positional argument, and [§3](@ref "3. One-sample") and
[§4](@ref "4. Two samples, equal variances (Student)") additionally accept summary
statistics (mean, standard deviation or variance, and count) in place of the data.

There is no departure from this specification to record.

## 10. References

Student's test is [Student (1908)](@cite student1908); the unequal-variance form and its
degrees of freedom are [Welch (1947)](@cite welch1947), following
[Satterthwaite (1946)](@cite satterthwaite1946).
