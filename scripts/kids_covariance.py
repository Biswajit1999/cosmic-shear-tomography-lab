"""Within-pair KiDS xi+ covariance and source-redshift diagnostics.

These are zero-shear null checks, not cosmological parameter fits. Pairwise
statistics must not be combined without the omitted cross-pair covariance.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import chi2


def _data(pair: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = np.asarray(pair['theta_arcmin'], dtype=float)
    x = np.asarray(pair['xi_plus'], dtype=float)
    c = np.asarray(pair['covariance'], dtype=float)
    if t.ndim != 1 or x.shape != t.shape or c.shape != (t.size, t.size):
        raise ValueError('incompatible angular, signal or covariance shapes')
    if t.size < 2 or np.any(t <= 0) or np.any(np.diff(t) <= 0):
        raise ValueError('angular bins must be increasing and positive')
    if not all(np.all(np.isfinite(z)) for z in (t, x, c)):
        raise ValueError('non-finite measurement or covariance')
    if not np.allclose(c, c.T, rtol=1e-10, atol=0):
        raise ValueError('non-symmetric covariance')
    return t, x, c


def _chi2(x: np.ndarray, c: np.ndarray) -> float:
    try:
        whitened = np.linalg.solve(np.linalg.cholesky(c), x)
    except np.linalg.LinAlgError as exc:
        raise ValueError('covariance must be positive definite') from exc
    return float(whitened @ whitened)


def eigenmodes(pair: dict) -> dict:
    """Null chi-square decomposed in correlation eigenmodes."""
    _, x, c = _data(pair)
    sigma = np.sqrt(np.diag(c))
    if np.any(sigma <= 0):
        raise ValueError('non-positive variance')
    eigenvalues, basis = np.linalg.eigh(c / np.outer(sigma, sigma))
    if np.any(eigenvalues <= 0):
        raise ValueError('non-positive covariance eigenvalue')
    contributions = (basis.T @ (x / sigma)) ** 2 / eigenvalues
    statistic = _chi2(x, c)
    if not np.isclose(sum(contributions), statistic, rtol=1e-8):
        raise ArithmeticError('eigenmodes disagree with Cholesky whitening')
    p = eigenvalues / eigenvalues.sum()
    return {'chi2': statistic, 'condition': float(eigenvalues[-1] / eigenvalues[0]),
            'effective_rank': float(np.exp(-sum(p * np.log(p)))),
            'mode_chi2': contributions.tolist()}


def angular_cut(pair: dict, minimum_arcmin: float, maximum_arcmin: float = np.inf) -> dict:
    """Marginal covariance submatrix for selected angular scales."""
    if not np.isfinite(minimum_arcmin) or np.isnan(maximum_arcmin) or minimum_arcmin > maximum_arcmin:
        raise ValueError('invalid angular cut bounds')
    theta, x, c = _data(pair)
    indices = np.flatnonzero((theta >= minimum_arcmin) & (theta <= maximum_arcmin))
    if not len(indices):
        raise ValueError('empty angular selection')
    stat = _chi2(x[indices], c[np.ix_(indices, indices)])
    return {'indices': indices.tolist(), 'chi2': stat, 'dof': int(len(indices)),
            'p_zero_shear': float(chi2.sf(stat, len(indices)))}


def nz_overlap(a: dict, b: dict) -> dict:
    """Bhattacharyya overlap of n(z), NOT lensing-kernel overlap."""
    def normal(item: dict) -> tuple[np.ndarray, np.ndarray]:
        z = np.asarray(item['z_mid'], dtype=float)
        n = np.asarray(item['density'], dtype=float)
        if z.ndim != 1 or n.shape != z.shape or len(z) < 2:
            raise ValueError('invalid n(z) shape')
        if np.any(~np.isfinite(z)) or np.any(~np.isfinite(n)) or np.any(n < 0):
            raise ValueError('invalid n(z) values')
        if np.any(np.diff(z) <= 0):
            raise ValueError('redshift grid must increase')
        integral = trapezoid(n, z)
        if integral <= 0:
            raise ValueError('zero n(z) integral')
        return z, n / integral

    za, na = normal(a)
    zb, nb = normal(b)
    grid = np.union1d(za, zb)
    pa = np.interp(grid, za, na, left=0, right=0)
    pb = np.interp(grid, zb, nb, left=0, right=0)
    pa /= trapezoid(pa, grid)
    pb /= trapezoid(pb, grid)
    overlap = float(np.clip(trapezoid(np.sqrt(pa * pb), grid), 0, 1))
    return {'bin_a': int(a['bin']), 'bin_b': int(b['bin']),
            'overlap': overlap, 'hellinger': float(np.sqrt(1 - overlap))}


def report(data: dict, minimum_arcmin: float = 6.0) -> dict:
    """Validate all released per-pair null statistics against the derivative."""
    if len(data['pairs']) != 15 or len(data['redshift_distributions']) != 5:
        raise ValueError('expected 15 tomographic pairs and five source bins')
    pairs = []
    for pair in data['pairs']:
        modes = eigenmodes(pair)
        if not np.isclose(modes['chi2'], pair['null_chi2'], rtol=1e-7):
            raise ValueError(f"inconsistent published chi2 for {pair['id']}")
        pairs.append({'id': pair['id'], 'eigen': modes,
                      'angular_cut': angular_cut(pair, minimum_arcmin)})
    bins = data['redshift_distributions']
    overlaps = [nz_overlap(a, b) for i, a in enumerate(bins) for b in bins[i:]]
    return {'source_fits_sha256': data['source_fits_sha256'],
            'cut_min_arcmin': minimum_arcmin, 'pairs': pairs, 'nz_overlaps': overlaps,
            'scope': 'Within-pair null statistics only; not a cosmological fit.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('data/kids1000-xip.json'))
    parser.add_argument('--output', type=Path)
    parser.add_argument('--min-arcmin', type=float, default=6.0)
    args = parser.parse_args()
    output = json.dumps(report(json.loads(args.input.read_text()), args.min_arcmin), indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output)
    else:
        print(output, end='')
