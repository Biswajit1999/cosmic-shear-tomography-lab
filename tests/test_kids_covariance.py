import numpy as np
import pytest
from scripts.kids_covariance import angular_cut, eigenmodes, nz_overlap, report


def pair():
    t = np.array([1, 2, 4, 8], float)
    s = np.array([1, 2, 3, 4], float) * 1e-4
    c = np.outer(s, s) * (0.2 + 0.8 * np.eye(4))
    x = np.array([1, -2, 0.5, 1]) * s
    return {'id': 'test', 'theta_arcmin': t.tolist(), 'xi_plus': x.tolist(),
            'covariance': c.tolist()}


def nz(b, mu):
    z = np.linspace(0, 3, 301)
    return {'bin': b, 'z_mid': z.tolist(),
            'density': np.exp(-0.5 * ((z - mu) / 0.2)**2).tolist()}


def test_eigenmodes_match_direct_chi2():
    p = pair()
    x, c = np.array(p['xi_plus']), np.array(p['covariance'])
    assert eigenmodes(p)['chi2'] == pytest.approx(x @ np.linalg.solve(c, x))
    assert sum(eigenmodes(p)['mode_chi2']) == pytest.approx(eigenmodes(p)['chi2'])


def test_scale_cut_uses_marginal_subcovariance():
    p = pair()
    x, c = np.array(p['xi_plus'])[2:], np.array(p['covariance'])[2:, 2:]
    assert angular_cut(p, 4)['chi2'] == pytest.approx(x @ np.linalg.solve(c, x))
    assert angular_cut(p, 0)['chi2'] == pytest.approx(eigenmodes(p)['chi2'])


def test_non_positive_covariance_rejected():
    p = pair()
    c = np.array(p['covariance'])
    c[0, 1] = c[1, 0] = 2 * np.sqrt(c[0, 0] * c[1, 1])
    p['covariance'] = c.tolist()
    with pytest.raises(ValueError):
        eigenmodes(p)


def test_redshift_overlap_limits_and_symmetry():
    a, b = nz(1, 0.5), nz(2, 2.5)
    assert nz_overlap(a, a)['overlap'] == pytest.approx(1)
    assert nz_overlap(a, b)['overlap'] < 0.001
    assert nz_overlap(a, b)['overlap'] == pytest.approx(nz_overlap(b, a)['overlap'])


def test_all_pairs_synthetic_report():
    pairs = []
    for i in range(15):
        p = pair()
        p['id'] = str(i)
        p['null_chi2'] = eigenmodes(p)['chi2']
        pairs.append(p)
    d = {'pairs': pairs, 'redshift_distributions': [nz(i, 0.3 * i) for i in range(1, 6)],
         'source_fits_sha256': 'synthetic-only'}
    result = report(d)
    assert len(result['pairs']) == 15
    assert len(result['nz_overlaps']) == 15
