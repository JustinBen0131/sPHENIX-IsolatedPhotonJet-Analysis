"""Validation: population, held-out performance and portable-model agreement.

The outer ten percent is evaluated after fitting, never used for epoch or
architecture selection. These numerical gates certify this model release,
not physics closure of a downstream analysis. Scores are evaluated in batches.
"""
from __future__ import annotations

import numpy as np

from .Models import require


def weighted_auc(y, score, weight):
    """Weighted rank AUC with tied scores kept together, as in the old trainer."""
    require(np.isfinite(score).all(), 'nonfinite prediction in validation')
    order = np.argsort(score, kind='stable')
    y, score, weight = y[order], score[order], weight[order]
    starts = np.r_[0, np.flatnonzero(np.diff(score)) + 1]
    positive = np.add.reduceat(weight * (y == 1), starts)
    negative = np.add.reduceat(weight * (y == 0), starts)
    require(positive.sum() > 0 and negative.sum() > 0, 'AUC needs both classes')
    previous = np.cumsum(negative) - negative
    return float((positive * (previous + .5 * negative)).sum() / (positive.sum() * negative.sum()))


def metrics(y, score, weight):
    """Weighted discrimination and probability losses for finite scores in [0, 1]."""
    require(np.isfinite(score).all() and ((score >= 0) & (score <= 1)).all(), 'score is outside [0,1]')
    clipped = np.clip(score, 1e-15, 1. - 1e-15)
    return dict(rows=len(y), classes=[int((y == k).sum()) for k in (0, 1)],
                auc=weighted_auc(y, score, weight),
                logloss=float(np.average(-y * np.log(clipped) - (1 - y) * np.log1p(-clipped), weights=weight)),
                brier=float(np.average((score - y) ** 2, weights=weight)))


def native_export(population, mask, reference, native, recipe, *, progress=None):
    """Check ROOT conversion against XGBoost on every applicable photon row.

    Float32 summation/sigmoid order can differ between the two libraries; the
    maximum absolute score error must be at most 1e-6. This export check does
    not change the tighter saved-native reload tolerance. Native scores own
    calibration and downstream decisions, avoiding mixed evaluator thresholds.
    """
    indices = np.flatnonzero(mask); maximum = 0.
    if progress: progress('native-export-validation', 0, len(indices), 'photons')
    tolerance = float(recipe['validation'].get('native_score_atol', 1e-6))
    require(0. < tolerance <= 1e-6, 'native score tolerance must be positive and at most 1e-6')
    for begin in range(0, len(indices), 65536):
        part = indices[begin:begin + 65536]
        original = np.asarray(reference(population.x[part]), np.float64)
        exported = np.asarray(native(population.x[part]), np.float64)
        require(original.shape == exported.shape == (len(part),), 'native export score shape differs')
        require(np.isfinite(original).all() and np.isfinite(exported).all() and
                ((original >= 0.) & (original <= 1.) & (exported >= 0.) & (exported <= 1.)).all(),
                'native export has invalid scores')
        maximum = max(maximum, float(np.max(np.abs(original - exported))) if len(part) else 0.)
        if progress: progress('native-export-validation', begin + len(part), len(indices), 'photons')
    require(len(indices) > 0 and maximum <= tolerance, 'native ROOT export differs from XGBoost')
    return dict(status='PASS', rows=len(indices), max_abs_error=maximum, absolute_tolerance=tolerance,
                feature_dtype='float32', calibration_evaluator='ROOT_RBDT')


def check(population, mask, trainer, saved, recipe, *, progress=None):
    """Compare every applicable saved score; report train and held-out metrics."""
    indices = np.flatnonzero(mask)
    score = np.empty(len(indices), np.float64); maximum = 0.
    if progress: progress('saved-score-validation', 0, len(indices), 'photons')
    for begin in range(0, len(indices), 65536):
        part = indices[begin:begin + 65536]
        original = np.asarray(trainer(population.x[part]), np.float64)
        reread = np.asarray(saved(population.x[part]), np.float64)
        require(original.shape == reread.shape == (len(part),), 'model score shape differs')
        require(np.isfinite(original).all() and np.isfinite(reread).all(), 'saved model has nonfinite scores')
        maximum = max(maximum, float(np.max(np.abs(original - reread))) if len(part) else 0.)
        score[begin:begin + len(part)] = reread
        if progress: progress('saved-score-validation', begin + len(part), len(indices), 'photons')
    require(maximum <= recipe['validation']['saved_score_atol'], 'saved model differs from trainer')
    train = population.train[indices]; holdout = population.holdout[indices]
    reports = {}
    for name, selected in [('train', train), ('holdout', holdout)]:
        if progress: progress('metrics-' + name, 0, 1, 'metric products')
        y, weight = population.y[indices[selected]], population.weights[indices[selected]]
        minimum = recipe['validation']['min_rows_per_class_per_split']
        require(all((y == k).sum() >= minimum for k in (0, 1)), name + ' lacks required class population')
        reports[name] = metrics(y, score[selected], weight)
        if progress: progress('metrics-' + name, 1, 1, 'metric products')
    if population.partition is not None and np.any(population.partition[indices] == 1):
        selected = population.partition[indices] == 1
        reports['validation'] = metrics(population.y[indices[selected]], score[selected], population.weights[indices[selected]])
    require(reports['holdout']['auc'] >= recipe['validation']['min_holdout_auc'], 'held-out discrimination failed')
    require(abs(reports['train']['auc'] - reports['holdout']['auc']) <=
            recipe['validation']['max_train_holdout_auc_gap'], 'train/holdout discrimination gap failed')
    return score, dict(status='PASS', saved_score_max_abs=maximum, split_sha256=population.receipt['split_sha256'],
                       tuning_on_holdout=False, **reports)
