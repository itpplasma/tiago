"""Independent observation-space checks of native grouped signal predictions."""
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.special import logsumexp, ndtr
from scipy.stats import multivariate_normal

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from grouped_signal_prediction import Noise, aggregate, log_predictive_density, predict


def main(executable):
    rng = np.random.default_rng(1702)
    x = np.c_[np.ones(12), rng.normal(size=(12, 2))]
    y = rng.normal(size=12)
    groups, times = np.repeat([0, 1, 2], 4), np.tile(np.arange(4)*.2, 3)
    new_x = np.c_[np.ones(6), rng.normal(size=(6, 2))]
    new_groups, new_times = np.repeat([8, 9], 3), np.tile(np.arange(3)*.3, 2)
    prior_sd = np.array([3., 1.7, .8])
    hypotheses = [Noise(.2), Noise(.3, .5, .8, .4), Noise(.5, .1, .4, 1.2)]
    prior_weights = np.array([.2, .5, .3])
    native = predict(executable, x, y, groups, times, new_x, new_groups, new_times,
                     prior_sd, hypotheses, prior_weights)
    oracle_means, oracle_covariances, logz = [], [], []
    for noise in hypotheses:
        def observation_covariance(a, ga, ta, b, gb, tb, same_rows):
            same_group = ga[:, None] == gb[None, :]
            result = (a*prior_sd**2)@b.T+noise.group_sd**2*same_group
            kernel = noise.residual_sd**2*noise.temporal_fraction*np.exp(
                -np.abs(ta[:, None]-tb[None, :])/noise.time_scale)*same_group
            result += kernel
            if same_rows:
                result += np.eye(len(a))*noise.residual_sd**2*(1-noise.temporal_fraction)
            return result
        covariance = observation_covariance(x, groups, times, x, groups, times, True)
        cross = observation_covariance(new_x, new_groups, new_times, x, groups, times, False)
        test = observation_covariance(new_x, new_groups, new_times,
                                      new_x, new_groups, new_times, True)
        factor = cho_factor(covariance)
        oracle_means.append(cross@cho_solve(factor, y))
        oracle_covariances.append(test-cross@cho_solve(factor, cross.T))
        logz.append(multivariate_normal.logpdf(y, cov=covariance))
    oracle_means, oracle_covariances = np.array(oracle_means), np.array(oracle_covariances)
    np.testing.assert_allclose(native['component_means'], oracle_means, atol=2e-10)
    np.testing.assert_allclose(native['component_covariances'], oracle_covariances, atol=2e-10)
    np.testing.assert_allclose(native['component_log_evidence'], logz, atol=2e-9)
    weights = np.exp(np.log(prior_weights)+logz-logsumexp(np.log(prior_weights)+logz))
    np.testing.assert_allclose(native['weights'], weights, atol=2e-10)
    for probability, bound in zip([.025, .975], native['interval95']):
        sd = np.sqrt(np.diagonal(oracle_covariances, axis1=1, axis2=2))
        np.testing.assert_allclose(weights@ndtr((bound-oracle_means)/sd), probability, atol=2e-10)
    values = rng.normal(size=6)
    expected_log = logsumexp([np.log(w)+multivariate_normal.logpdf(values, m, c)
                             for w, m, c in zip(weights, oracle_means, oracle_covariances)])
    np.testing.assert_allclose(log_predictive_density(native, values), expected_log, atol=2e-9)
    average = aggregate(native, [0, 1, 2])
    np.testing.assert_allclose(average['component_variances'],
                              oracle_covariances[:, :3, :3].sum(axis=(1, 2))/9, atol=2e-10)
    # A new group's common offset cannot average away with more samples.
    many = predict(executable, x, y, groups, times, np.zeros((30, 3)),
                   np.full(30, 7), np.arange(30), prior_sd, [Noise(.1, .6)])
    variance = aggregate(many, np.arange(30))['component_variances'][0]
    np.testing.assert_allclose(variance, .6**2+.1**2/30, atol=2e-12)
    for bad_groups in [groups[:6], np.repeat([0, 8], 3)]:
        try:
            predict(executable, x, y, groups, times, new_x, bad_groups, new_times,
                    prior_sd, hypotheses)
        except ValueError:
            pass
        else:
            raise AssertionError('overlapping prediction groups accepted')
    for labels in [np.full(6, np.nan), np.full(6, np.datetime64('NaT', 'ns'))]:
        try:
            predict(executable, x, y, groups, times, new_x, labels, new_times,
                    prior_sd, hypotheses)
        except ValueError:
            pass
        else:
            raise AssertionError('missing group labels accepted')
    try:
        aggregate(native, [0.9])
    except ValueError:
        pass
    else:
        raise AssertionError('fractional aggregation indices accepted')
    for parameters in [(0., 0., 0., .5), (.1, -.1, 0., .5), (.1, .1, 1., .5)]:
        try:
            Noise(*parameters)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid covariance prior accepted')
    print('Observation-space conditioning/evidence, mixture quantiles, joint density and group averages: PASS')


if __name__ == '__main__':
    main(sys.argv[1])
