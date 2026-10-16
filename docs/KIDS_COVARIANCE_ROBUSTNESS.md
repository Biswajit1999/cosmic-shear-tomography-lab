# KiDS-1000 xi+ covariance robustness and angular-shape diagnostics

Run `python scripts/analyze_shear_robustness.py` against the **already committed** `data/kids1000-xip.json`. The script produces `results/kids_covariance_robustness.json` without network access. Run `python -m unittest discover -s tests -p 'test_shear_diagnostics.py' -v` from the repository root.

Three numerical diagnostic families are implemented:

1. **Covariance audit:** Cholesky-whiten each nine-angle xi+ vector using its released **within-pair** 9×9 covariance. Report χ²_full, χ²_diagonal-only, the correlation-matrix condition number, participation-ratio effective rank, and largest off-diagonal correlation. The diagonal approximation is a sensitivity comparator, not the official likelihood.
2. **Leave-one-angle-out conditional leverage:** For each angular bin, recompute χ² on the eight-angle submatrix and evaluate Δχ². The Schur complement guarantees Δχ² = (conditional residual / conditional sigma)² ≥ 0 for a positive-definite covariance. The nine values are correlated; do not treat them as nine independent significances or jackknife uncertainty estimates.
3. **Covariance-weighted angular profile:** Grid-search a phenomenological xi+(θ)=A(θ/10 arcmin)^(-γ) profile, solving for A analytically at each γ. A is unconstrained and can be negative. This is a descriptive angular-shape comparison, not a physically motivated matter-power-spectrum model, S8 measurement, or calibrated detection statistic.

**Critical limitation:** The browser derivative preserves each pair's 9×9 covariance but omits the cross-pair covariance among all 15 tomographic pairs. The output must never be used to combine pair-level p-values as independent measurements or to claim a KiDS cosmological posterior. No existing FITS source, browser data, or published null test is changed.
