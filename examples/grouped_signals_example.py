"""Independent-shot synthetic check of shared-offset predictive uncertainty."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from grouped_signal_prediction import Noise, aggregate, noise_prior, predict


def main(executable, output):
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261107)
    count, train_groups, test_groups = 8, 10, 80
    group = np.repeat(np.arange(train_groups+test_groups), count)
    times = np.tile(np.arange(count)*.5, train_groups+test_groups)
    features = rng.normal(size=(len(group), 2))
    x = np.c_[np.ones(len(group)), features]
    beta = np.array([1.5, .7, -.3])
    truth_noise = Noise(.08, .35, .75, .5)
    offset = rng.normal(scale=truth_noise.group_sd, size=train_groups+test_groups)
    residual = np.empty(len(group))
    for label in range(train_groups+test_groups):
        take = group == label
        covariance = truth_noise.covariance(group[take], times[take])
        residual[take] = np.linalg.cholesky(covariance)@rng.normal(size=count)
    y = x@beta+offset[group]+residual
    train, test = group < train_groups, group >= train_groups
    hypotheses, masses = noise_prior(refinement=1)
    records, arrays = {}, {}
    for model in ['no_shot_offset', 'shot_offset_mixture']:
        take = np.array([h.group_sd == 0 or model == 'shot_offset_mixture' for h in hypotheses])
        result = predict(executable, x[train], y[train], group[train], times[train],
            x[test], group[test], times[test], np.array([5., 2., 2.]),
            [h for h, keep in zip(hypotheses, take) if keep], masses[take])
        target = y[test]
        interval = result['interval95']
        averages = [aggregate(result, np.flatnonzero(group[test] == label))
                    for label in range(train_groups, train_groups+test_groups)]
        recorded_means = np.array([target[group[test] == label].mean()
                                   for label in range(train_groups, train_groups+test_groups)])
        predicted_means = np.array([a['mean'] for a in averages])
        mean_intervals = np.array([a['interval95'] for a in averages]).T
        records[model] = dict(pointwise_rmse=float(np.sqrt(np.mean((target-result['mean'])**2))),
            covered_blocks=int(np.sum((target >= interval[0]) & (target <= interval[1]))),
            total_blocks=len(target), covered_shot_means=int(np.sum(
                (recorded_means >= mean_intervals[0]) & (recorded_means <= mean_intervals[1]))),
            total_shots=test_groups, mean_interval_width=float(np.mean(mean_intervals[1]-mean_intervals[0])),
            block_interval_width=float(np.mean(interval[1]-interval[0])),
            shot_mean_rmse=float(np.sqrt(np.mean((recorded_means-predicted_means)**2))))
        arrays[model+'_recorded'] = recorded_means
        arrays[model+'_predicted'] = predicted_means
        arrays[model+'_interval95'] = mean_intervals
    records['provenance'] = dict(seed=20261107, training_groups=train_groups, test_groups=test_groups,
        blocks_per_group=count, coefficients=beta.tolist(),
        truth=dict(residual_sd=.08, group_sd=.35, temporal_fraction=.75, time_scale=.5),
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        qualification='Independent synthetic groups under a fixed true relation; empirical coverage is not a general coverage guarantee')
    (output/'results.json').write_text(json.dumps(records, indent=2)+'\n')
    np.savez_compressed(output/'predictions.npz', **arrays)
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.executable.resolve(), args.output)
