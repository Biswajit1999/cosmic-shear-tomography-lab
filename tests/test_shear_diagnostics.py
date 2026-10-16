"""Synthetic matrix invariants and exact analytic special cases."""
import unittest
import numpy as np
from research_analysis.shear_diagnostics import (
    covariance_audit, leave_one_angle_out, power_law_profile, validated_pair)


def fixture(rho=0.0):
    theta = np.array([1, 2, 4, 8, 16.], float)
    y = np.array([.3, .1, -.2, .4, -.1])
    corr = (1-rho)*np.eye(5) + rho*np.ones((5, 5))
    sig = np.array([.2, .3, .4, .5, .6])
    cov = corr*np.outer(sig, sig)
    return {'theta_arcmin': theta.tolist(), 'xi_plus': y.tolist(),
            'covariance': cov.tolist()}


class CovarianceTests(unittest.TestCase):
    def test_diagonal_equals_full(self):
        a = covariance_audit(fixture())
        self.assertAlmostEqual(a['chi2_full'], a['chi2_diagonal_only'])
        self.assertAlmostEqual(a['correlation_effective_rank'], 5)

    def test_nonzero_correlations(self):
        a = covariance_audit(fixture(.4))
        self.assertGreater(a['max_abs_offdiag_correlation'], .39)
        self.assertGreater(a['correlation_condition_number'], 1)
        self.assertNotAlmostEqual(a['chi2_full'], a['chi2_diagonal_only'])

    def test_amplitude_scaling(self):
        p = fixture(.3)
        a = covariance_audit(p)['chi2_full']
        p['xi_plus'] = (np.asarray(p['xi_plus'])*2).tolist()
        b = covariance_audit(p)['chi2_full']
        self.assertAlmostEqual(b, 4*a)

    def test_covariance_scaling(self):
        p = fixture(.3)
        a = covariance_audit(p)['chi2_full']
        p['covariance'] = (2*np.asarray(p['covariance'])).tolist()
        self.assertAlmostEqual(covariance_audit(p)['chi2_full'], a/2)

    def test_reject_non_spd(self):
        p = fixture()
        p['covariance'][0][0] = -1
        with self.assertRaises(ValueError): validated_pair(p)

    def test_reject_non_monotonic_theta(self):
        p = fixture()
        p['theta_arcmin'][1] = 0.5
        with self.assertRaises(ValueError): validated_pair(p)


class ConditionalTests(unittest.TestCase):
    def test_diagonal_exact_contributions(self):
        p = fixture()
        y = np.array(p['xi_plus'])
        var = np.diag(p['covariance'])
        result = leave_one_angle_out(p)
        for e, expectation in zip(result['angles'], y*y/var):
            self.assertAlmostEqual(e['conditional_delta_chi2'], expectation)

    def test_correlated_schur_identity(self):
        p = fixture(.65)
        result = leave_one_angle_out(p)
        for e in result['angles']:
            self.assertGreaterEqual(e['conditional_variance'], 0)
            self.assertAlmostEqual(e['conditional_delta_chi2'], e['conditional_z']**2, places=8)

    def test_zero_signal(self):
        p = fixture(.5)
        p['xi_plus'] = [0]*5
        self.assertEqual(leave_one_angle_out(p)['max_conditional_delta_chi2'], 0)


class PowerLawTests(unittest.TestCase):
    def test_exact_power_law(self):
        p = fixture(.25)
        theta = np.asarray(p['theta_arcmin'])
        p['xi_plus'] = (3e-5*(theta/10)**(-1.2)).tolist()
        result = power_law_profile(p, gamma_grid=np.linspace(.5, 1.8, 131))
        self.assertAlmostEqual(result['best']['gamma'], 1.2, places=10)
        self.assertAlmostEqual(result['best']['amplitude'], 3e-5, places=10)
        self.assertLess(result['best']['chi2_residual'], 1e-22)

    def test_amplitude_sign(self):
        p = fixture()
        theta = np.asarray(p['theta_arcmin'])
        p['xi_plus'] = (-1e-5*(theta/10)**(-.8)).tolist()
        result = power_law_profile(p, gamma_grid=np.linspace(.5, 1, 51))
        self.assertLess(result['best']['amplitude'], 0)

    def test_reject_bad_grid(self):
        with self.assertRaises(ValueError): power_law_profile(fixture(), gamma_grid=np.array([1., 0.]))

if __name__ == '__main__': unittest.main()
