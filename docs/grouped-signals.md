# Predicting signals across shots

Repeated samples within one shot share conditions and possible calibration
offsets. Treating them as independent shots can make the uncertainty of a shot
average too small. `tools/grouped_signal_prediction.py` composes TIAGO's native
linear Gaussian inference with group offsets, correlated residuals and a finite
Bayesian prior over their strengths. It predicts recorded signals in wholly
unseen groups, rather than noiseless plasma quantities.

## Model

For observation `i`, the signal is

\[
y_i=x_i^T\beta+b_{g(i)}+\epsilon_i.
\]

Here `x_i` is the recorded predictor vector, `beta` the shared regression
coefficients, and `g(i)` its shot or acquisition-group label. Coefficients have
independent zero-mean Gaussian priors with user-specified positive standard
deviations. Include a column of ones for an intercept. Each group offset `b`
is an independent zero-mean Gaussian with standard deviation `tau`.

The residual `epsilon` has covariance

\[
R_{ij}=\sigma^2[(1-\rho)\delta_{ij}
  +\rho\,1_{g(i)=g(j)}\exp(-|t_i-t_j|/\ell)].
\]

`sigma > 0` is the residual standard deviation, `0 <= rho < 1` its temporal
fraction, `t` the acquisition time and `ell > 0` the correlation time in the
same units. The Kronecker delta `delta` is one on the diagonal and zero
elsewhere. Different groups have independent residuals and offsets. This is
a declared covariance family, not a guarantee that real drift follows it.

Each `Noise` hypothesis fixes these four quantities. Its regression coefficients
and training-group offsets are integrated analytically by the native
`tiago_linear_infer` executable, retaining likelihood normalization. Posterior
hypothesis masses are proportional to prior mass times marginal likelihood.
Predictions form a Gaussian mixture, including coefficient uncertainty, a new
group's offset, residual covariance and uncertainty between hypotheses.
These are exact Gaussian conditionals, not approximations to nonlinear modes.

The group offset does not average away as more points from the same group are
combined. Prediction covariance also couples different new groups through
their uncertain shared coefficients. Quantiles come from the mixture cumulative
distribution, rather than a Gaussian approximation to its moments.

## Python interface

Add the checkout's `tools` directory to the Python path. NumPy and SciPy are
required, together with a built `tiago_linear_infer` executable.

```python
from grouped_signal_prediction import Noise, aggregate, predict

result = predict(executable, x, y, groups, times, new_x, new_groups, new_times,
                 prior_sd, [Noise(.1), Noise(.1, .3, .75, .5)], [.25, .75])
average = aggregate(result, [0, 1, 2])
```

Training designs have observations in rows, and `prior_sd` matches their
columns. Group labels are integers or nonempty strings. Every prediction group
must be absent from training. Within-group prediction would require residual
cross-covariance and conditioning on that group's existing observations; this
interface rejects it. Predictor values are conditioned on as known. Their
measurement error, missing predictor entries and calibration need a larger
generative model.

Results include component masses, means, full covariance matrices, normalized
training evidences, mixture moments and pointwise 95% predictive intervals.
`aggregate` computes a selected group's equal-weight mean distribution using
full covariance. `log_predictive_density(result, values, indices=None)` gives
the joint predictive log density, optionally marginalized to selected indices.
Group-mean intervals are distinct from simultaneous bands for all points.
Component covariance storage grows quadratically with prediction count; request
modest batches when full covariance across many groups is unnecessary.

`noise_prior` constructs normalized finite priors using log-amplitude trapezoid
weights, a separate zero-offset atom and two temporal fractions. Each resulting
finite mixture is a proper discrete Bayesian model. Increasing `refinement`
approximates integration under the stated continuous log-uniform amplitudes;
check convergence before interpreting tiny evidence or predictive differences.
Changing prior support is a separate scientific sensitivity check.

## Validation and interpretation

`tests/check_grouped_prediction.py` compares native coefficient-space inference
with independent observation-space Gaussian conditioning and marginal densities.
It checks mixture quantiles, joint prediction, temporal covariance, unseen-group
aggregation and a common offset that cannot vanish through averaging. Missing
group identifiers and fractional aggregation indices are rejected.

`examples/grouped_signals_example.py` trains on independent synthetic shots and
withholds entire new shots. Its shared-offset experiment compares pointwise
intervals with intervals for shot means. This example prescribes offsets
between shots; it establishes no particular diagnostic's short-time drift.

The conditional Gaussian machinery is standard Bayesian regression, described
in [Rasmussen and Williams, chapter 2](https://gaussianprocess.org/gpml/chapters/RW2.pdf).
Validate actual transfers across configurations and dates. A correctly computed
posterior can fail when that model does not describe the withheld experiment.
More coverage obtained merely by widening intervals is not sufficient: inspect
prediction error, interval width and properly correlated predictive density.

Missing signals remain posterior predictions. Never add them as independent
measurements alongside their predictors, or multiply this conditional likelihood
by another likelihood for the same target data. Combining learned relationships
with physical inference requires one explicit joint generative model and
qualified shared calibration/processing. This module alone does not identify
magnetic islands or certify experimental uncertainty.
