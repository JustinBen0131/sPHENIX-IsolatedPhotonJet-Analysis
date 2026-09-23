"""Small numerical predecessors, copied verbatim. See HISTORICAL_SOURCES.json.

Only pure quantile/fit functions are retained here; no historical orchestration.
The isolation HistogramData is the original record required by its quantile.
"""
from dataclasses import dataclass
import math
import numpy as np

def weighted_threshold(scores: np.ndarray, weights: np.ndarray, efficiency: float) -> float:
    """Verbatim THE-134 sealed quantile: stable mergesort, ties never split, exact sample value."""
    scores = np.asarray(scores, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    good = np.isfinite(scores) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(good):
        return math.nan
    scores = scores[good]
    weights = weights[good]
    order = np.argsort(scores, kind="mergesort")
    scores = scores[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    target = (1.0 - efficiency) * float(cumulative[-1])
    index = int(np.searchsorted(cumulative, target, side="left"))
    return float(scores[min(index, len(scores) - 1)])

def weighted_efficiency(mask: np.ndarray, weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=np.float64)
    denominator = float(np.sum(weights))
    return float(np.sum(weights[np.asarray(mask, dtype=bool)]) / denominator) if denominator > 0 else math.nan

def threshold_uncertainty(scores: np.ndarray, weights: np.ndarray, threshold: float,
                          efficiency: float, *, density_window: float,
                          score_histogram_bins: int) -> dict[str, float]:
    """THE-96 quantile uncertainty transcribed onto the score axis.

    Kish effective entries give the binomial error on the quantile probability, and dividing by the
    local weighted density converts it to an error on the threshold. The density is the symmetric
    finite difference of the weighted quantile function over ``density_window`` in probability, which
    is the stable form of the single-bin density THE-96 reads off its isolation histogram.
    """
    total = float(np.sum(weights))
    sumw2 = float(np.sum(np.square(weights)))
    width = 1.0 / float(score_histogram_bins)
    effective_entries = total * total / sumw2 if sumw2 > 0.0 else total
    lower_efficiency = min(efficiency + density_window, 1.0)
    upper_efficiency = max(efficiency - density_window, 0.0)
    lower = weighted_threshold(scores, weights, lower_efficiency)
    upper = weighted_threshold(scores, weights, upper_efficiency)
    span = float(upper - lower)
    probability_span = float(lower_efficiency - upper_efficiency)
    density = probability_span / span if span > 0.0 else 0.0
    if effective_entries > 0.0 and density > 0.0:
        probability_error = math.sqrt(efficiency * (1.0 - efficiency) / effective_entries)
        statistical_error = probability_error / density
    else:
        statistical_error = width
    return {
        "threshold_error_gev": float(max(0.5 * width, statistical_error)),
        "sumw": total,
        "sumw2": sumw2,
        "effective_entries": float(effective_entries),
        "local_score_density": float(density),
    }

def constant_fit(y: np.ndarray, yerr: np.ndarray) -> dict[str, float]:
    """Verbatim THE-96 weighted constant fit."""
    valid = np.isfinite(y) & np.isfinite(yerr) & (yerr > 0.0)
    if int(np.count_nonzero(valid)) < 2:
        raise ValueError("Constant fit needs at least two valid points")
    y = y[valid]
    yerr = yerr[valid]
    weights = 1.0 / np.square(yerr)
    value = float(np.sum(weights * y) / np.sum(weights))
    error = float(math.sqrt(1.0 / np.sum(weights)))
    chi2 = float(np.sum(np.square((y - value) / yerr)))
    return {"value_gev": value, "error_gev": error, "chi2": chi2, "ndf": int(y.size - 1), "n_points": int(y.size)}

def linear_fit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray) -> dict:
    """Verbatim THE-96 weighted linear fit."""
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0.0)
    if int(np.count_nonzero(valid)) < 3:
        raise ValueError("Linear fit needs at least three valid points")
    x = x[valid]
    y = y[valid]
    yerr = yerr[valid]
    design = np.column_stack([np.ones_like(x), x])
    weights = 1.0 / np.square(yerr)
    normal = design.T @ (weights[:, None] * design)
    covariance = np.linalg.inv(normal)
    parameters = covariance @ (design.T @ (weights * y))
    prediction = design @ parameters
    chi2 = float(np.sum(np.square((y - prediction) / yerr)))
    return {
        "intercept_gev": float(parameters[0]),
        "slope_gev_per_percent": float(parameters[1]),
        "intercept_error_gev": float(math.sqrt(covariance[0, 0])),
        "slope_error_gev_per_percent": float(math.sqrt(covariance[1, 1])),
        "covariance": covariance.tolist(),
        "chi2": chi2,
        "ndf": int(y.size - 2),
        "n_points": int(y.size),
    }

@dataclass(frozen=True)
class HistogramData:
    values: np.ndarray
    variances: np.ndarray
    edges: np.ndarray
    underflow: float
    overflow: float
    underflow_variance: float
    overflow_variance: float
    sources: tuple[str, ...]

def quantile_from_histogram(data: HistogramData, probability: float) -> dict[str, float]:
    values = np.asarray(data.values, dtype=float)
    negative_tolerance = 1e-10 * max(1.0, float(np.max(np.abs(values))))
    if np.any(values < -negative_tolerance):
        raise ValueError("Isolation histogram has negative weighted bins; CDF quantile is undefined")
    values = np.maximum(values, 0.0)
    in_range_total = float(np.sum(values))
    underflow = max(float(data.underflow), 0.0)
    overflow = max(float(data.overflow), 0.0)
    total = in_range_total + underflow + overflow
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("Isolation histogram has no positive in-range weight")

    cumulative = np.cumsum(values)
    absolute_target = probability * total
    if absolute_target <= underflow:
        raise ValueError(
            f"Requested quantile {probability:.3f} lies in histogram underflow below {data.edges[0]:g} GeV"
        )
    if absolute_target > underflow + in_range_total:
        raise ValueError(
            f"Requested quantile {probability:.3f} lies in histogram overflow above {data.edges[-1]:g} GeV"
        )
    target = absolute_target - underflow
    index = int(np.searchsorted(cumulative, target, side="left"))
    index = min(max(index, 0), len(values) - 1)
    previous = float(cumulative[index - 1]) if index > 0 else 0.0
    bin_weight = float(values[index])
    fraction = 0.5 if bin_weight <= 0.0 else (target - previous) / bin_weight
    fraction = min(max(float(fraction), 0.0), 1.0)
    width = float(data.edges[index + 1] - data.edges[index])
    threshold = float(data.edges[index] + fraction * width)

    sumw2 = float(
        np.sum(np.maximum(data.variances, 0.0))
        + data.underflow_variance
        + data.overflow_variance
    )
    effective_entries = total * total / sumw2 if sumw2 > 0.0 else total
    density = bin_weight / (total * width) if width > 0.0 else 0.0
    if effective_entries > 0.0 and density > 0.0:
        probability_error = math.sqrt(probability * (1.0 - probability) / effective_entries)
        statistical_error = probability_error / density
    else:
        statistical_error = width
    threshold_error = max(0.5 * width, statistical_error)
    flow_fraction = (underflow + overflow) / max(total, 1e-12)
    return {
        "threshold_gev": threshold,
        "threshold_error_gev": float(threshold_error),
        "sumw": total,
        "sumw2": sumw2,
        "effective_entries": float(effective_entries),
        "underflow": data.underflow,
        "overflow": data.overflow,
        "flow_fraction": float(flow_fraction),
    }
