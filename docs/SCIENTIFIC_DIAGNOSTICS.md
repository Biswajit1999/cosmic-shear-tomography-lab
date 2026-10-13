# KiDS-1000 scientific covariance diagnostics

This analysis uses the committed KiDS-1000 xi+ derivative and adds three independent numerical checks.

**Correlation eigenmodes:** whiten the nine-bin covariance, diagonalize its correlation matrix, and verify that the sum of mode contributions reproduces the Cholesky-based zero-shear chi-square. Report the correlation condition number and entropy effective rank.

**Angular scale cuts:** choose bins by angular separation and use the corresponding *marginal covariance submatrix* to calculate the zero-shear chi-square. A subblock of the inverse full covariance is not equivalent and must not be substituted.

**Source n(z) overlap:** calculate normalized Bhattacharyya overlap and Hellinger distance for the five released redshift distributions. These quantify redshift distribution similarity, not lensing-kernel overlap.

Run `python -m pytest tests/test_kids_covariance.py -q` for synthetic numerical tests and `python scripts/kids_covariance.py --input data/kids1000-xip.json --output results/kids_covariance.json` for reproducible survey-derivative checks.

The analysis preserves the FITS source hash and validates all 15 stored within-pair chi-square values. It does not combine pairwise p-values, fit cosmological parameters, infer S8, or model the omitted cross-pair covariance.
