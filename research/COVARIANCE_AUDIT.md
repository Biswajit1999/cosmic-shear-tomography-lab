# KiDS within-pair covariance sensitivity (22 November 2026)

Two independent conditional calculations: (1) compare the zero-shear quadratic form as the released within-pair covariance is interpolated toward its diagonal; (2) profile a purely phenomenological constant xi+ offset with generalized least squares. Input: one `pairs[]` entry from `data/kids1000-xip.json` (`xi_plus`, `covariance`). Run from `research`: `python -m unittest -v test_covariance_audit.py`.

The interpolation is a robustness scenario, not an estimated covariance correction. A constant offset is not a shear theory. The JSON excludes cross-tomographic-pair covariance: do not combine pairwise significances or claim S8 constraints. Tests use synthetic numbers, not observational inference. Code and tests prepared with AI assistance for review under Biswajit Jana's scientific responsibility.
