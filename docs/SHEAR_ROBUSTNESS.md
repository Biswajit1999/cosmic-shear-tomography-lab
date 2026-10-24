# KiDS-1000 within-pair robustness diagnostics

Run from the repository root after installing `requirements-research.txt`:

```bash
python -m pytest -q tests/test_shear_robustness.py
python scripts/run_shear_robustness.py --output shear_robustness_results.json
```

This independent Python diagnostic reads **only** the existing `data/kids1000-xip.json` derivative, which records the KiDS-1000 `xi+` values, nine angular bins, and their 9×9 within-pair covariance for each of 15 tomographic pairs. The source FITS hash is carried into the report. No observations are fabricated or replaced.

Three checks address distinct numerical questions: (1) correlation-matrix eigenvalues, conditioning, and the effect of dropping off-diagonal covariance; (2) leave-one-angular-bin-out zero-shear chi-square using each exact principal covariance submatrix, whose difference is a squared *conditional innovation* and is not additive across bins; (3) an empirical power-law angular profile with analytically profiled amplitude and a finite slope grid. The power-law model is descriptive, **not** a physical cosmic-shear theory prediction; its nominal degrees of freedom are not a validated goodness-of-fit likelihood.

**Limitations:** the browser derivative does not retain cross-tomographic-pair covariance. Consequently, these statistics must not be combined into a joint 135-bin detection significance, S8 posterior, or cosmological parameter constraint. The covariance is treated as supplied and fixed; uncertainty in its estimation is not modelled. Angle exclusion changes the statistical question, so comparing p-values after cuts is a sensitivity study rather than an independent detection test.

Data citation: Asgari et al. (2021), *A&A* 645 A104; KiDS-1000 public cosmic-shear two-point-correlation release, as linked in the repository README.
