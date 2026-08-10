# Cosmic Shear Tomography Lab

Weak-lensing correlation functions and redshift-bin sensitivity.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of the simulation. The app loads `data/reference.json`, renders those published anchors first, then sends the adjustable model to `physicsWorker.js` so numerical work stays off the UI thread.

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation.
- `data/reference.json`: small auditable reference-data bundle.
- `scripts/validate.js`: no-dependency repository validation.

## Run

```bash
python -m http.server 8080
```

Open `http://localhost:8080`.

## Validate

```bash
npm run check
```

The validation script checks required files, JSON reference data, worker syntax, citations and absence of unfinished scaffold tokens.

## Reference Data

Compressed xi_plus-like weak-lensing anchors for validating tomographic trend behaviour.

## References

- Bartelmann, M. and Schneider, P., 2001. Weak gravitational lensing. Physics Reports, 340(4-5), pp.291-472.
- Troxel, M.A. et al., 2018. Dark Energy Survey Year 1 results: Cosmological constraints from cosmic shear. Physical Review D, 98(4), p.043528.

## Research Quality Upgrade

See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for the validation layer, reference anchors, equations and research boundaries added to this repository.
