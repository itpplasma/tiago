# Predicting an unrecorded signal

Use a small joint Gaussian model when simultaneous training observations support
a roughly linear relationship between recorded signals. A proper
normal-inverse-Wishart prior represents uncertainty in their common mean and
covariance. Integrating both parameters gives a multivariate Student-t
posterior predictive distribution. Conditioning on the available channels gives
a linear prediction for the missing channels, with uncertainty from finite
training data. No sampling or neural network is required.

This model predicts recorded signals. It does not by itself separate plasma
variation, measurement noise, calibration drift or model discrepancy. Those
components require an explicit observation model and independent information.
Use training-only transformations and retain their units and parameters.

## Fortran interface

`tiago_signal_correlation::signal_correlation_t` provides:

```fortran
call model%fit(samples, mean0, psi0, kappa0, nu0, status)
call model%predict(observed, values, missing, location, scale, df, log_observed, status)
```

`samples(d,n)` contains `n` complete training vectors of dimension `d`.
`mean0(d)`, `psi0(d,d)`, `kappa0 > 0` and `nu0 > d-1` specify the prior.
`psi0` must be symmetric positive definite. Its convention is
`E[covariance] = psi0 / (nu0-d-1)` when that expectation exists. Zero training
vectors are supported. Finite values are required; zero status means success.

`observed` contains distinct one-based channel indices in any order, with
matching `values`. Outputs list missing indices in ascending order, their
conditional location and Student-t scale matrix, degrees of freedom `df`, and
the normalized log marginal density of the observed values. The scale matrix
is **not the covariance**: the latter is `df/(df-2)*scale` for `df > 2` and
does not exist otherwise. The conditional distribution remains defined.
Observing no channels returns the joint predictive distribution and zero log
marginal density. Observing every channel returns empty prediction arrays.

## Batch executable

`tiago_signal_predict INPUT OUTPUT` fits once and evaluates each query
separately with OpenMP. The input contains:

1. `d ntraining nqueries kappa0 nu0`.
2. One prior mean row, followed by `d` rows of `psi0`.
3. One row per complete training vector.
4. For each query: the observed count, then an index row and a value row.
   Omit the last two rows when the count is zero.

For each query, output begins with
`query_index missing_count df log_observed`. If channels are missing, this is
followed by their indices, conditional locations and one scale row per missing
channel. `examples/signal_prediction.py` supplies a Python caller and a
reproducible withheld-sensor demonstration:

```sh
python3 examples/signal_prediction.py --executable build/tiago_signal_predict --output /tmp/signals
```

The example needs NumPy and SciPy. Its stable-sensor experiment is accompanied
by changes in the withheld sensor's offset and correlation. A good prediction
under stable conditions does not guarantee transfer under those changes.

## Use inside inference

A predicted channel is a posterior predictive quantity, not another independent
measurement. Adding it to the likelihood alongside the channels used to predict
it counts their information twice. Use the joint likelihood marginalized to the
actually observed channels instead. This class supplies that observed marginal
for its own joint-signal model; it does not automatically combine it with a
separate physical likelihood that uses the same observations.

If independently labelled classes have different joint-signal models, their
observed marginal densities can be compared with explicit class priors.
If a class changes only a missing channel and the observed marginal is identical,
the available signals cannot update the class probabilities. An unconditional
correlation model alone does not identify magnetic topology.

Validate on entire withheld shots or independent acquisition blocks. Treating
many correlated samples as independent training vectors makes parameter
uncertainty too small. Challenge altered calibration, configuration and
operating conditions, and record whether missingness depends on those
conditions. Report predictive coverage and error against a training-only
marginal baseline. The present model does not model temporal dependence or
partially observed training rows; neither is silently imputed.

The conjugate update follows
[Murphy's Gaussian conjugate analysis](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf).
Conditional Student distributions are also described in
[Shah, Wilson and Ghahramani (2014)](https://arxiv.org/abs/1402.4306).
Tests compare the batch update with sequential updates and conditional densities
with an independent joint-density/marginal-density ratio.
