"""WorkingPoints: weighted full-signal efficiency thresholds for each model.

Reuse the existing weighted signal quantiles and linear threshold convention.
AuAu fits versus centrality; optional pp fits versus photon ET. Routed models
use their own cell calibration, never another model's thresholds. One-bin
centrality cells have a constant local threshold. Strict score comparisons
retain ties and require ordered, non-overlapping ID regions.
"""
from __future__ import annotations

import numpy as np

from .Models import require


def weighted_quantile(values, weights, q):
    """Return an observed score at the weighted CDF crossing, without splitting ties."""
    require(len(values) > 0 and np.isfinite(values).all() and np.isfinite(weights).all() and
            (weights > 0).all(), 'invalid working-point signal population')
    order = np.argsort(values, kind='stable')
    values, weights = values[order], weights[order]
    cdf = np.cumsum(weights)
    index = int(np.searchsorted(cdf, max(0., min(1., q)) * cdf[-1], side='left'))
    return float(values[min(index, len(values) - 1)])


def derive(population, mask, score, cell, recipe):
    """Calibrate WP70/80/90 on every eligible signal row in this cell.

    Train and holdout both enter calibration. Return the centrality/ET line
    and its measured efficiencies; insufficient bins or failed strict-cut
    efficiency/order checks refuse the release.
    """
    indices = np.flatnonzero(mask)
    signal = population.y[indices] == 1
    values = population.centrality[indices] if recipe['system'] == 'AuAu' else population.et[indices]
    weights = population.weights[indices]
    config = recipe['working_points']
    axis = 'cent_value' if recipe['system'] == 'AuAu' else 'pho_et'
    edges = config['centrality_edges' if recipe['system'] == 'AuAu' else 'et_edges']
    bounds = [0., 80.] if recipe['system'] == 'AuAu' else [15., 35.]
    for key, limits in cell['bounds'].items():
        if (recipe['system'] == 'AuAu' and key.startswith('centrality')) or (recipe['system'] == 'PP' and key == 'et'):
            bounds = limits
    edges = sorted(set([bounds[0], bounds[1]] + [e for e in edges if bounds[0] < e < bounds[1]]))
    working_points = {}
    for efficiency in config['efficiencies']:
        rows = []
        for low, high in zip(edges[:-1], edges[1:]):
            selected = signal & (values >= low) & (values < high)
            require(int(selected.sum()) >= 20, 'working-point calibration bin has too few signal photons')
            threshold = weighted_quantile(score[selected], weights[selected], 1. - efficiency)
            exact_efficiency = float(np.average(score[selected] > threshold, weights=weights[selected]))
            require(abs(exact_efficiency - efficiency) <= config['max_efficiency_error'],
                    'strict weighted-quantile efficiency gate failed; score ties cannot be split')
            rows.append(dict(low=low, high=high, center=(low + high) / 2.,
                             threshold=threshold, signal=int(selected.sum()), exact_efficiency=exact_efficiency))
        centers = np.array([r['center'] for r in rows]); thresholds = np.array([r['threshold'] for r in rows])
        if len(rows) == 1:
            intercept, slope = float(thresholds[0]), 0.
        else:
            # The maintained working-point representation is an affine line.
            slope, intercept = np.polyfit(centers, thresholds, 1)
            intercept, slope = float(intercept), float(slope)
        residual = thresholds - (intercept + slope * centers)
        errors = []
        for row in rows:
            selected = signal & (values >= row['low']) & (values < row['high'])
            passed = score[selected] > intercept + slope * values[selected]
            achieved = float(np.average(passed, weights=weights[selected]))
            row['achieved_efficiency'] = achieved
            errors.append(abs(achieved - efficiency))
        require(max(errors) <= config['max_efficiency_error'], 'working-point efficiency gate failed')
        require(float(np.max(np.abs(residual))) <= config['max_fit_residual'], 'working-point fit residual failed')
        require(float(np.sqrt(np.mean(residual ** 2))) <= config['max_fit_rms'], 'working-point fit RMS failed')
        working_points['WP' + str(round(efficiency * 100))] = dict(
            intercept=intercept, slope=slope, target=efficiency, bins=rows,
            max_efficiency_error=max(errors), fit_rms=float(np.sqrt(np.mean(residual ** 2))))
    require(set(working_points) == {'WP70', 'WP80', 'WP90'}, 'required working points missing')
    for point in bounds:
        cuts = {k: v['intercept'] + v['slope'] * point for k, v in working_points.items()}
        require(.1 < cuts['WP80'] <= cuts['WP70'] <= 1., 'tight/non-tight regions overlap or are empty')
        require(0. <= cuts['WP90'] <= cuts['WP80'], 'working-point efficiency order differs')
    return dict(axis=axis, population='all eligible signal, train plus holdout',
                signal_photons=int(signal.sum()), bounds=bounds, cuts=working_points, status='PASS')


def classify(score, axis_value, working_points):
    """Return ID codes: -1 unscored, 0 below, 1 non-tight, 2 gap, 3 tight."""
    cuts = working_points['cuts']
    tight = cuts['WP70']['intercept'] + cuts['WP70']['slope'] * axis_value
    high = cuts['WP80']['intercept'] + cuts['WP80']['slope'] * axis_value
    return np.where(np.isfinite(score), np.where(score > tight, 3,
                    np.where((score > .1) & (score < high), 1, np.where(score >= high, 2, 0))), -1)
