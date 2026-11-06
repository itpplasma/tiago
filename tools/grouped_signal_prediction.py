"""Bayesian linear signal prediction with group offsets and time correlation."""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import brentq
from scipy.special import logsumexp, ndtr, ndtri
from scipy.stats import multivariate_normal

from linear_batch import infer


def noise_prior(residual_bounds=(.025, .8), group_bounds=(.05, .8),
                refinement=0, zero_mass=.25, time_scale=.5):
    """Normalized finite priors; log-amplitude trapezoids plus a zero-offset atom."""
    if (not isinstance(refinement, int) or refinement < 0 or not 0 < zero_mass < 1):
        raise ValueError('invalid noise-prior refinement or atom')
    nodes, masses = [], []
    for bounds, intervals in [(residual_bounds, 5), (group_bounds, 4)]:
        if len(bounds) != 2 or not 0 < bounds[0] < bounds[1] or not np.all(np.isfinite(bounds)):
            raise ValueError('invalid positive amplitude bounds')
        count = intervals*2**refinement+1
        nodes.append(np.geomspace(*bounds, count))
        weights = np.ones(count)
        weights[[0, -1]] = .5
        masses.append(weights/np.sum(weights))
    groups = np.r_[0., nodes[1]]
    group_masses = np.r_[zero_mass, (1-zero_mass)*masses[1]]
    hypotheses, weights = [], []
    for sd, sd_mass in zip(nodes[0], masses[0]):
        for tau, tau_mass in zip(groups, group_masses):
            for fraction in [0., .75]:
                hypotheses.append(Noise(sd, tau, fraction, time_scale))
                weights.append(sd_mass*tau_mass/2)
    return hypotheses, np.array(weights)


def _group_labels(values):
    labels = np.asarray(values)
    if labels.dtype.kind not in 'iuUS' or (labels.dtype.kind in 'US' and np.any(labels == '')):
        raise ValueError('group labels must be integers or nonempty strings')
    return labels


@dataclass(frozen=True)
class Noise:
    residual_sd: float
    group_sd: float = 0.
    temporal_fraction: float = 0.
    time_scale: float = .5

    def __post_init__(self):
        values = [self.residual_sd, self.group_sd, self.temporal_fraction, self.time_scale]
        if (not np.all(np.isfinite(values)) or self.residual_sd <= 0
                or self.group_sd < 0 or not 0 <= self.temporal_fraction < 1
                or self.time_scale <= 0):
            raise ValueError('invalid group/time noise parameters')

    def covariance(self, groups, times):
        groups = _group_labels(groups)
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or groups.shape != times.shape or not np.all(np.isfinite(times)):
            raise ValueError("invalid group/time coordinates")
        same = groups[:, None] == groups[None, :]
        temporal = np.exp(-abs(times[:, None]-times[None, :])/self.time_scale)*same
        return self.residual_sd**2*((1-self.temporal_fraction)*np.eye(len(times))
                                   + self.temporal_fraction*temporal)


def mixture_quantiles(means, variances, weights, probabilities=(.025, .975)):
    """Pointwise quantiles; variances have component and prediction axes."""
    sd = np.sqrt(variances)
    result = np.empty((len(probabilities), means.shape[1]))
    for row, probability in enumerate(probabilities):
        if not 0 < probability < 1:
            raise ValueError('invalid quantile probability')
        bounds = means+ndtri(probability)*sd
        for column in range(means.shape[1]):
            lo, hi = np.min(bounds[:, column]), np.max(bounds[:, column])
            def residual(value):
                return np.dot(weights, ndtr((value-means[:, column])/sd[:, column]))-probability
            result[row, column] = lo if lo == hi else brentq(residual, lo, hi)
    return result


def predict(executable, x, y, groups, times, new_x, new_groups, new_times,
            prior_sd, hypotheses, prior_weights=None):
    """Predict complete new groups; all continuous fit coordinates are integrated."""
    x, y, times, new_x, new_times, prior_sd = [np.asarray(v, dtype=float) for v in
        [x, y, times, new_x, new_times, prior_sd]]
    groups, new_groups = _group_labels(groups), _group_labels(new_groups)
    if x.ndim != 2 or new_x.ndim != 2 or x.shape[1] != new_x.shape[1]:
        raise ValueError('training/prediction design mismatch')
    n, p = x.shape
    nn = len(new_x)
    if (n < 1 or nn < 1 or p < 1 or y.shape != (n,) or times.shape != (n,)
            or groups.shape != (n,) or new_times.shape != (nn,)
            or new_groups.shape != (nn,) or prior_sd.shape != (p,)):
        raise ValueError('invalid grouped prediction dimensions')
    if any(not np.all(np.isfinite(v)) for v in [x, y, times, new_x, new_times, prior_sd]):
        raise ValueError('nonfinite grouped prediction input')
    if np.any(prior_sd <= 0) or np.intersect1d(groups, new_groups).size:
        raise ValueError('positive coefficient priors and wholly new groups required')
    if not hypotheses:
        raise ValueError('empty noise prior')
    weights0 = (np.ones(len(hypotheses)) if prior_weights is None
                else np.asarray(prior_weights, dtype=float))
    if (weights0.shape != (len(hypotheses),) or not np.all(np.isfinite(weights0))
            or np.any(weights0 <= 0)):
        raise ValueError('invalid hypothesis prior masses')
    log_prior = np.log(weights0)-logsumexp(np.log(weights0))
    labels = np.unique(groups)
    membership = (groups[:, None] == labels).astype(float)
    batch = []
    for noise in hypotheses:
        design = np.c_[x, membership] if noise.group_sd else x
        sd = np.r_[prior_sd, np.full(len(labels), noise.group_sd)] if noise.group_sd else prior_sd
        lower = cholesky(noise.covariance(groups, times), lower=True)
        batch.append(dict(design=solve_triangular(lower, design, lower=True),
                          data=solve_triangular(lower, y, lower=True), prior_sd=sd,
                          log_normalization=-n/2*np.log(2*np.pi)-np.log(np.diag(lower)).sum()))
    fits = infer(executable, batch)
    log_evidence = np.array([fit['log_evidence'] for fit in fits])
    log_total = logsumexp(log_prior+log_evidence)
    weights = np.exp(log_prior+log_evidence-log_total)
    means, covariances = [], []
    same_new = new_groups[:, None] == new_groups[None, :]
    for noise, fit in zip(hypotheses, fits):
        means.append(new_x@fit['mean'][:p])
        covariance = new_x@fit['covariance'][:p, :p]@new_x.T
        covariance += noise.group_sd**2*same_new+noise.covariance(new_groups, new_times)
        covariances.append((covariance+covariance.T)/2)
    means, covariances = np.array(means), np.array(covariances)
    mean = weights@means
    centred = means-mean
    covariance = np.einsum('h,hij->ij', weights, covariances)
    covariance += np.einsum('h,hi,hj->ij', weights, centred, centred)
    intervals = mixture_quantiles(means, np.diagonal(covariances, axis1=1, axis2=2), weights)
    return dict(mean=mean, covariance=covariance, interval95=intervals, weights=weights,
                component_means=means, component_covariances=covariances,
                component_log_evidence=log_evidence, log_evidence=float(log_total))


def aggregate(result, indices):
    """Distribution of the equally weighted mean over selected predictions."""
    indices = np.asarray(indices)
    if indices.dtype.kind not in 'iu':
        raise ValueError('aggregation indices must be integers')
    if (indices.ndim != 1 or not len(indices) or np.any(indices < 0)
            or np.any(indices >= result['mean'].size) or len(np.unique(indices)) != len(indices)):
        raise ValueError('invalid aggregation indices')
    means = result['component_means'][:, indices].mean(axis=1)[:, None]
    cov = result['component_covariances'][:, indices][:, :, indices]
    variances = cov.mean(axis=(1, 2))[:, None]
    interval = mixture_quantiles(means, variances, result['weights'])[:, 0]
    return dict(mean=float(result['weights']@means[:, 0]), interval95=interval,
                component_means=means[:, 0], component_variances=variances[:, 0])


def log_predictive_density(result, values, indices=None):
    """Joint predictive density; includes temporal and shared-group covariance."""
    means, covariances = result['component_means'], result['component_covariances']
    if indices is not None:
        indices = np.asarray(indices)
        if indices.ndim != 1 or indices.dtype.kind not in 'iu' or not len(indices):
            raise ValueError('marginal indices must be a nonempty integer vector')
        if np.any(indices < 0) or np.any(indices >= result['mean'].size) or len(np.unique(indices)) != len(indices):
            raise ValueError('invalid marginal indices')
        means, covariances = means[:, indices], covariances[:, indices][:, :, indices]
    return float(logsumexp([np.log(weight)+multivariate_normal.logpdf(values, mean, covariance)
        for weight, mean, covariance in zip(result['weights'], means, covariances) if weight > 0]))
