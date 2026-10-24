"""Analytical and physics-aware regression tests for within-pair diagnostics."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from shear_robustness import covariance_diagnostic, leave_one_angle_out, null_diagnostic, powerlaw_profile


def pair(y=None, rho=0):
    theta = np.geomspace(1, 100, 9)
    sigma = 0.02 * np.ones(9)
    corr = rho ** np.abs(np.subtract.outer(np.arange(9), np.arange(9)))
    return {"theta_arcmin": theta, "xi_plus": np.asarray(y if y is not None else np.ones(9) * 0.04),
            "covariance": corr * np.outer(sigma, sigma)}


def test_diagonal_null_exact():
    p = pair()
    assert null_diagnostic(p)["chi2"] == pytest.approx(36)
    assert covariance_diagnostic(p)["diagonal_only_chi2"] == pytest.approx(36)


def test_correlated_null_matches_inverse():
    p = pair(rho=0.65)
    y = np.asarray(p["xi_plus"])
    assert null_diagnostic(p)["chi2"] == pytest.approx(y @ np.linalg.solve(p["covariance"], y))


def test_angular_principal_submatrix_not_diagonal():
    p = pair(rho=0.8)
    indices = [0, 3, 6]
    c = p["covariance"][np.ix_(indices, indices)]
    y = p["xi_plus"][indices]
    assert null_diagnostic(p, indices)["chi2"] == pytest.approx(y @ np.linalg.solve(c, y))


def test_leave_one_out_conditional_identity():
    p = pair(y=np.arange(1, 10) * 0.01, rho=0.4)
    full = null_diagnostic(p)["chi2"]
    loo = leave_one_angle_out(p)
    for i, result in enumerate(loo):
        reduced = null_diagnostic(p, [j for j in range(9) if j != i])["chi2"]
        assert result["delta_chi2"] == pytest.approx(full - reduced, abs=1e-9)
        assert result["delta_chi2"] >= 0


def test_powerlaw_recovers_slope_and_amplitude():
    p = pair(rho=0.35)
    theta = p["theta_arcmin"]
    pivot = np.exp(np.mean(np.log(theta)))
    p["xi_plus"] = 0.06 * (theta / pivot) ** (-1.25)
    result = powerlaw_profile(p, np.linspace(0, 2.5, 501))
    assert result["alpha"] == pytest.approx(1.25, abs=0.005)
    assert result["amplitude"] == pytest.approx(0.06, rel=1e-8)
    assert result["chi2"] < 1e-17


def test_powerlaw_amplitude_standard_error():
    p = pair()
    result = powerlaw_profile(p, [0, 0.5])
    assert result["amplitude_sigma_fixed_alpha"] > 0


def test_singular_covariance_rejected():
    p = pair()
    p["covariance"] = np.ones((9, 9))
    with pytest.raises(ValueError, match="positive definite"):
        null_diagnostic(p)


def test_nonsymmetric_covariance_rejected():
    p = pair()
    p["covariance"][0, 1] = 1e-5
    with pytest.raises(ValueError, match="symmetric"):
        covariance_diagnostic(p)


def test_invalid_angular_index_rejected():
    with pytest.raises(ValueError, match="out of range"):
        null_diagnostic(pair(), [9])


def test_zero_signal_null_is_zero():
    p = pair(y=np.zeros(9), rho=0.9)
    assert null_diagnostic(p)["chi2"] == pytest.approx(0)
    assert all(x["delta_chi2"] == pytest.approx(0) for x in leave_one_angle_out(p))
