"""Within-pair KiDS-1000 xi+ diagnostics with explicit covariance limitations.

No inter-pair covariance is stored in the browser derivative. Never combine
pair-level p-values or interpret these diagnostics as a cosmological likelihood.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.stats import chi2


def validated_pair(pair: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Check finite xi+, strictly ordered positive angles and SPD covariance."""
    y = np.asarray(pair['xi_plus'], dtype=float)
    theta = np.asarray(pair['theta_arcmin'], dtype=float)
    cov = np.asarray(pair['covariance'], dtype=float)
    if y.ndim != 1 or y.size < 3 or theta.shape != y.shape or cov.shape != (y.size, y.size):
        raise ValueError('inconsistent vector/angle/covariance dimensions')
    if not np.isfinite(y).all() or not np.isfinite(theta).all() or not np.isfinite(cov).all():
        raise ValueError('nonfinite values')
    if np.any(theta <= 0) or np.any(np.diff(theta) <= 0):
        raise ValueError('angles must be strictly increasing and positive')
    if not np.allclose(cov, cov.T, atol=0, rtol=1e-10):
        raise ValueError('covariance must be symmetric')
    try:
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as e:
        raise ValueError('covariance must be positive definite') from e
    return theta, y, cov


def covariance_audit(pair: dict) -> dict:
    """Eigen-spectrum, correlation structure, and full/diagonal null tests."""
    theta, y, cov = validated_pair(pair)
    sigma = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sigma, sigma)
    eigen = np.linalg.eigvalsh(corr)
    whitened = np.linalg.solve(np.linalg.cholesky(cov), y)
    full_chi2 = float(whitened @ whitened)
    diagonal_chi2 = float(np.sum((y / sigma)**2))
    offdiag = corr - np.diag(np.diag(corr))
    return {
        'n_angles': len(theta),
        'chi2_full': full_chi2,
        'chi2_diagonal_only': diagonal_chi2,
        'p_full_null': float(chi2.sf(full_chi2, len(theta))),
        'correlated_distance': math.sqrt(full_chi2),
        'correlation_condition_number': float(eigen[-1]/eigen[0]),
        'correlation_effective_rank': float(eigen.sum()**2 / np.dot(eigen, eigen)),
        'max_abs_offdiag_correlation': float(np.max(np.abs(offdiag))),
        'whitened_residuals': whitened.tolist(),
        'warning': 'single tomographic pair only; diagonal-only statistic is not the released likelihood',
    }


def leave_one_angle_out(pair: dict) -> dict:
    """Conditional contribution of each angle to the correlated null chi-square.

    Delta = chi2_all - chi2_without_i >= 0 by the Schur complement.
    These are *not* independent significances or jackknife error bars.
    """
    theta, y, cov = validated_pair(pair)
    total = float(y @ np.linalg.solve(cov, y))
    entries = []
    for i, angle in enumerate(theta):
        keep = np.arange(y.size) != i
        y0 = y[keep]
        c0 = cov[np.ix_(keep, keep)]
        removed_chi2 = float(y0 @ np.linalg.solve(c0, y0))
        delta = max(0.0, total - removed_chi2)  # clip roundoff only
        # Conditional residual and variance from the Schur complement.
        c_cross = cov[i, keep]
        predicted = float(c_cross @ np.linalg.solve(c0, y0))
        conditional_var = float(cov[i, i] - c_cross @ np.linalg.solve(c0, c_cross))
        standardized = (float(y[i]) - predicted) / math.sqrt(conditional_var)
        entries.append({
            'theta_arcmin': float(angle), 'chi2_without_angle': removed_chi2,
            'conditional_delta_chi2': delta,
            'conditional_z': standardized,
            'conditional_variance': conditional_var,
        })
    return {'chi2_all': total, 'angles': entries,
            'max_conditional_delta_chi2': max(e['conditional_delta_chi2'] for e in entries),
            'warning': 'conditional leverage diagnostics, not independent detections'}


def power_law_profile(pair: dict, gamma_grid: np.ndarray | None = None,
                      pivot_arcmin: float = 10.0) -> dict:
    """Covariance-weighted A(theta/pivot)^(-gamma) descriptive shape fit.

    A is unconstrained and can be negative. Grid search over gamma introduces
    model-selection effects: returned chi2 must not be treated as a cosmology fit.
    """
    theta, y, cov = validated_pair(pair)
    if not math.isfinite(pivot_arcmin) or pivot_arcmin <= 0:
        raise ValueError('pivot must be positive')
    gamma = np.linspace(-0.5, 2.5, 121) if gamma_grid is None else np.asarray(gamma_grid, dtype=float)
    if gamma.ndim != 1 or gamma.size < 2 or not np.isfinite(gamma).all() or np.any(np.diff(gamma) <= 0):
        raise ValueError('gamma grid must be finite, increasing, and length >= 2')
    L = np.linalg.cholesky(cov)
    yw = np.linalg.solve(L, y)
    profiles = []
    for g in gamma:
        template = (theta / pivot_arcmin) ** (-g)
        tw = np.linalg.solve(L, template)
        norm = float(tw @ tw)
        amplitude = float(tw @ yw / norm)
        sigma = 1 / math.sqrt(norm)
        resid = yw - amplitude * tw
        profiles.append({'gamma': float(g), 'amplitude': amplitude,
                         'amplitude_sigma': sigma,
                         'chi2_residual': float(resid @ resid)})
    best = min(profiles, key=lambda x: x['chi2_residual'])
    return {'pivot_arcmin': float(pivot_arcmin), 'best': best,
            'n_angles': len(theta), 'nominal_residual_dof': len(theta)-2,
            'profile': profiles,
            'warning': 'phenomenological angular-shape diagnostic, not S8 inference'}
