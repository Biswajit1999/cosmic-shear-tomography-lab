"""Covariance-aware diagnostics for the released KiDS-1000 xi+ browser derivative.

These are per-tomographic-pair checks, not an S8 or cosmology likelihood.
Only the 9x9 within-pair covariance retained in data/kids1000-xip.json is used.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import chi2


def _pair_arrays(pair):
    theta = np.asarray(pair["theta_arcmin"], dtype=float)
    signal = np.asarray(pair["xi_plus"], dtype=float)
    cov = np.asarray(pair["covariance"], dtype=float)
    n = len(theta)
    if n < 3 or signal.shape != (n,) or cov.shape != (n, n):
        raise ValueError("inconsistent pair dimensions or too few angular bins")
    if not (np.all(np.isfinite(theta)) and np.all(np.isfinite(signal)) and np.all(np.isfinite(cov))):
        raise ValueError("all inputs must be finite")
    if np.any(theta <= 0) or np.any(np.diff(theta) <= 0):
        raise ValueError("angles must be positive and strictly increasing")
    if not np.allclose(cov, cov.T, rtol=1e-11, atol=0):
        raise ValueError("covariance must be symmetric")
    if np.any(np.diag(cov) <= 0):
        raise ValueError("covariance diagonal must be positive")
    scale = np.sqrt(np.diag(cov))
    corr = cov / np.outer(scale, scale)
    try:
        np.linalg.cholesky(corr)
    except np.linalg.LinAlgError as exc:
        raise ValueError("covariance must be positive definite") from exc
    return theta, signal, cov


def _quadratic(signal, covariance):
    scale = np.sqrt(np.diag(covariance))
    corr = covariance / np.outer(scale, scale)
    z = signal / scale
    chol = np.linalg.cholesky(corr)
    return float(np.dot(np.linalg.solve(chol, z), np.linalg.solve(chol, z)))


def null_diagnostic(pair, indices=None):
    """Test xi+=0 using the selected principal covariance submatrix."""
    theta, signal, cov = _pair_arrays(pair)
    idx = np.arange(len(theta)) if indices is None else np.asarray(indices, dtype=int)
    if idx.ndim != 1 or len(idx) < 1 or len(set(idx.tolist())) != len(idx):
        raise ValueError("indices must be a nonempty unique 1-D selection")
    if np.any(idx < 0) or np.any(idx >= len(theta)):
        raise ValueError("angular-bin index out of range")
    c = cov[np.ix_(idx, idx)]
    chi_square = _quadratic(signal[idx], c)
    return {"chi2": chi_square, "dof": int(len(idx)),
            "p_zero_shear": float(chi2.sf(chi_square, len(idx))),
            "correlated_distance": float(np.sqrt(chi_square))}


def covariance_diagnostic(pair):
    """Report dimensionless correlation eigenvalues and Cholesky whitening."""
    _, y, cov = _pair_arrays(pair)
    sigma = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sigma, sigma)
    eigenvalues = np.linalg.eigvalsh(corr)
    return {"correlation_condition": float(eigenvalues[-1] / eigenvalues[0]),
            "min_correlation_eigenvalue": float(eigenvalues[0]),
            "max_correlation_eigenvalue": float(eigenvalues[-1]),
            "null_chi2": _quadratic(y, cov),
            "diagonal_only_chi2": float(np.sum((y / sigma) ** 2))}


def leave_one_angle_out(pair):
    """Conditional chi2 contributions: full minus reduced-principal chi2.

    The difference is a squared conditional innovation, not an independent
    detection statistic. Individual contributions are not additive.
    """
    theta, _, _ = _pair_arrays(pair)
    full = null_diagnostic(pair)["chi2"]
    result = []
    for i, angle in enumerate(theta):
        reduced = null_diagnostic(pair, np.delete(np.arange(len(theta)), i))["chi2"]
        delta = full - reduced
        if delta < -1e-8 * max(1.0, full):
            raise ArithmeticError("conditional chi-square contribution is negative")
        result.append({"angle_arcmin": float(angle), "delta_chi2": float(max(0, delta)),
                       "reduced_chi2": reduced})
    return result


def powerlaw_profile(pair, alpha_grid=None):
    """Profile a descriptive xi+(theta)=A*(theta/pivot)^(-alpha) template.

    This empirical template is not a theoretical weak-lensing prediction.
    It uses the complete within-pair covariance and profiles amplitude exactly.
    """
    theta, signal, cov = _pair_arrays(pair)
    alpha = np.linspace(0, 2.5, 251) if alpha_grid is None else np.asarray(alpha_grid, dtype=float)
    if alpha.ndim != 1 or len(alpha) < 2 or not np.all(np.isfinite(alpha)):
        raise ValueError("alpha grid must contain at least two finite values")
    if np.any(np.diff(alpha) <= 0) or np.any(alpha < 0) or np.any(alpha > 4):
        raise ValueError("alpha grid must be increasing within [0,4]")
    pivot = float(np.exp(np.mean(np.log(theta))))
    template = (theta[None, :] / pivot) ** (-alpha[:, None])
    scale = np.sqrt(np.diag(cov))
    chol = np.linalg.cholesky(cov / np.outer(scale, scale))
    whitened_t = np.linalg.solve(chol, (template / scale).T)
    whitened_y = np.linalg.solve(chol, signal / scale)
    numerator = whitened_t.T @ whitened_y
    denominator = np.sum(whitened_t**2, axis=0)
    amplitudes = numerator / denominator
    chi2_grid = np.sum((whitened_y[:, None] - whitened_t * amplitudes[None, :])**2, axis=0)
    index = int(np.argmin(chi2_grid))
    return {"alpha": float(alpha[index]), "amplitude": float(amplitudes[index]),
            "amplitude_sigma_fixed_alpha": float(1 / np.sqrt(denominator[index])),
            "chi2": float(chi2_grid[index]), "dof_nominal": int(len(theta) - 2),
            "pivot_arcmin": pivot, "grid_boundary_minimum": bool(index in (0, len(alpha) - 1))}
