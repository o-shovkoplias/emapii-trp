"""Block bootstrap and Van't Hoff regression."""
from __future__ import annotations

from collections.abc import Callable

import numpy as np

from . import R_KJ


def block_bootstrap(
    x: np.ndarray,
    stat: Callable[[np.ndarray], float],
    *,
    block: int,
    n_boot: int,
    rng: np.random.Generator,
    ci: float = 0.95,
) -> tuple[float, float, float]:
    """Moving-block bootstrap of ``stat`` over a time series ``x``.

    Returns (point estimate, lower CI bound, upper CI bound). Blocks of ``block``
    consecutive frames are resampled with replacement until the original length is
    reached, which preserves the short-time autocorrelation of the series.
    """
    x = np.asarray(x)
    n = len(x)
    block = max(1, min(block, n))
    n_blocks = int(np.ceil(n / block))
    starts_max = n - block + 1
    vals = np.empty(n_boot)
    for b in range(n_boot):
        starts = rng.integers(0, starts_max, n_blocks)
        idx = (starts[:, None] + np.arange(block)[None, :]).ravel()[:n]
        vals[b] = stat(x[idx])
    a = (1 - ci) / 2
    return float(stat(x)), float(np.nanquantile(vals, a)), float(np.nanquantile(vals, 1 - a))


def integrated_autocorr_time(x: np.ndarray, max_lag: int | None = None) -> float:
    """Integrated autocorrelation time (frames) of a 0/1 indicator by summing the ACF
    up to its first zero crossing (Sokal-style initial-positive-sequence estimator)."""
    x = np.asarray(x, float)
    x = x - x.mean()
    if x.var() == 0:
        return float("nan")
    n = len(x)
    max_lag = max_lag or n // 4
    acf = np.correlate(x, x, mode="full")[n - 1:n - 1 + max_lag] / (x.var() * n)
    tau = 1.0
    for k in range(1, max_lag):
        if acf[k] <= 0:
            break
        tau += 2 * acf[k]
    return float(tau)


def vant_hoff_fit(T_K: np.ndarray, lnK: np.ndarray) -> dict:
    """Least-squares fit ln K = -dH/(R T) + dS/R.

    Returns dH (kJ/mol), dS (J/mol/K), T_c = dH/dS (K and C), R^2, n, and the standard
    error of the slope/intercept (assumes independent points; n <= 6 here so treat as
    indicative only).
    """
    T_K = np.asarray(T_K, float)
    lnK = np.asarray(lnK, float)
    x = 1.0 / T_K
    n = len(x)
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, lnK, rcond=None)
    fit = slope * x + intercept
    ss_res = float(np.sum((lnK - fit) ** 2))
    ss_tot = float(np.sum((lnK - lnK.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    dH = -slope * R_KJ
    dS = intercept * R_KJ
    Tc_K = dH / dS if abs(dS) > 1e-12 else float("nan")
    if n > 2:
        sigma2 = ss_res / (n - 2)
        cov = np.linalg.inv(A.T @ A) * sigma2
        se_slope, se_int = np.sqrt(np.diag(cov))
    else:
        se_slope = se_int = float("nan")
    return {
        "n": int(n), "slope": float(slope), "intercept": float(intercept),
        "dH_kJmol": float(dH), "dS_JmolK": float(dS * 1000), "Tc_K": float(Tc_K),
        "Tc_C": float(Tc_K - 273.15), "R2": float(r2),
        "se_dH_kJmol": float(se_slope * R_KJ), "se_dS_JmolK": float(se_int * R_KJ * 1000),
    }
