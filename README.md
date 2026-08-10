# Cosmic Shear Tomography Lab

Weak-lensing correlation functions and redshift-bin sensitivity, explored interactively in the browser.

Created and maintained by Biswajit Jana.

## Scientific Background

### Weak gravitational lensing

Light from distant galaxies is deflected by the intervening distribution of matter (dark and
baryonic) as it travels to us. This deflection subtly distorts the observed shapes of background
("source") galaxies — stretching them tangentially around foreground mass concentrations. In the
weak regime (away from cluster cores and strong-lensing arcs) this distortion is a few percent
effect, far smaller than the intrinsic scatter in any single galaxy's shape. It is only measurable
statistically, by correlating the shapes of many galaxy pairs across the sky.

### Cosmic shear as a large-scale-structure probe

"Cosmic shear" is the name for weak lensing produced not by individual halos but by the full,
continuous distribution of large-scale structure between us and a broad population of source
galaxies. Because the lensing signal depends directly on the projected matter distribution, cosmic
shear is one of the cleanest probes of the amplitude and growth of structure, largely independent
of galaxy bias. It is a pillar statistic for Stage-III/IV surveys (KiDS, DES, HSC, and the
upcoming Euclid and LSST/Rubin surveys), used to constrain the matter density Ω_m and the
clustering amplitude σ8.

### Shear correlation functions ξ+ / ξ−

The two-point shear correlation functions ξ+(θ) and ξ−(θ) are the standard real-space statistics
used to summarize the cosmic shear field as a function of angular separation θ between galaxy
pairs. They are constructed from the tangential and cross ellipticity components of galaxy pairs
and are related to the convergence power spectrum P_κ(ℓ) through Bessel-function projections
(see the Math Appendix). ξ+ is the sum of the tangential and cross correlations and is sensitive
to a broad range of angular scales; ξ− isolates smaller scales and is more sensitive to
non-linear structure growth. This lab focuses on a ξ+-like observable.

### Tomographic binning

Because every source galaxy has some (photometric) redshift estimate, source galaxies can be split
into tomographic bins by estimated redshift. Correlating shear within and across bins isolates the
lensing efficiency and growth-of-structure information at different cosmic epochs, dramatically
increasing the constraining power of a cosmic shear survey relative to a single, redshift-averaged
measurement. The mean source redshift of a bin, z_source, sets the lensing kernel: higher-redshift
bins probe lensing from mass along a longer path length and typically show a stronger shear signal.

### Connection to S8 and σ8

Cosmic shear amplitude scales approximately as σ8·√(Ω_m / 0.3), a combination commonly reported as
S8 = σ8√(Ω_m/0.3). Comparing the S8 value inferred from cosmic shear surveys against the value
inferred from Planck CMB observations has been a long-running tension in the field (the "S8
tension"), making the amplitude scaling implemented in this lab directly relevant to current
cosmology debates.

### Survey context

Real-world cosmic shear measurements come from wide-field imaging surveys: the Kilo-Degree Survey
(KiDS), the Dark Energy Survey (DES), and — imminently — Euclid and the Vera C. Rubin
Observatory's LSST. These surveys measure ξ+/ξ− (or equivalent band-power/COSEBI statistics)
across tomographic bins and fit them against theoretical predictions marginalized over intrinsic
alignments, baryonic feedback, and photometric-redshift uncertainty. The reference anchor points
bundled with this lab (`data/reference.json`) are illustrative points in that same observable
space (θ, ξ+), used only to give the interactive model a fixed visual and numerical anchor — they
are not a reproduction of any survey's actual measured data vector.

## How It Works

This is a zero-build, static browser laboratory — there is no backend or build step.

1. `index.html` lays out a "mission control" dashboard: a controls panel, a combined
   reference-data/simulation plot, a parameter-response heatmap, and a telemetry readout.
2. On load, `app.js` fetches `data/reference.json` and plots the bundled reference anchor points
   (θ, ξ+) directly on the canvas, alongside their citation and dataset description.
3. `app.js` then posts the current model parameters to `physicsWorker.js`, a Web Worker, so the
   numerical evaluation (curve sampling + heatmap generation) never blocks the UI thread.
4. `physicsWorker.js` evaluates a closed-form ξ+-like model (see the `shear()` function) as a
   function of angular separation θ, using four adjustable parameters:
   - `sigma8` — the matter clustering amplitude σ8
   - `omegaM` — the matter density parameter Ω_m
   - `zSource` — the mean redshift of the tomographic source bin
   - `ia` — an intrinsic-alignment amplitude that suppresses the observed signal
   The worker returns a sampled ξ+(θ) curve, a set of summary metrics (ξ+ at θ=10 arcmin, the
   derived S8 = σ8√(Ω_m/0.3), and the IA suppression factor), and a 2D "tomographic kernel"
   heatmap used purely as a qualitative visualization of how the signal amplitude responds to σ8
   and bin depth.
5. Moving any control slider re-runs the worker in real time; the plot, heatmap, and telemetry
   panel update accordingly, and the reference anchor points stay fixed on the plot as a visual
   check on the model's scale and shape.
6. `research-overlay.js` adds a small, non-invasive quality panel showing validation status and
   the benchmark anchors in `data/research-reference.json`, which are checked (independent of the
   physics model) by `scripts/validate_repository.mjs`.

This is a teaching/exploration tool, not a cosmological inference pipeline: the model is a
compact analytic approximation chosen to reproduce the qualitative scaling behaviour of ξ+ (power-law
decline at small scales, exponential large-scale suppression, σ8/Ω_m/redshift/IA dependence), not
a full Limber-integrated theory prediction fit to real survey data.

## Usage

```bash
python -m http.server 8080
```

Open `http://localhost:8080` and adjust the sliders (σ8, Ω_m, source-bin mean z, intrinsic-alignment
amplitude) to see the ξ+(θ) curve, the tomographic-kernel heatmap, and the derived S8/telemetry
values respond in real time. Use "Reset model" to return to the fiducial parameter values.

## Validate

```bash
npm run check
```

The validation script checks required files, JSON reference data, worker syntax, citations, and
absence of unfinished scaffold tokens. See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for details
on the validation layer, reference anchors, and research boundaries.

## Architecture

- `index.html` — mission-control interface.
- `styles.css` — dense dark scientific dashboard.
- `app.js` — UI state, Canvas rendering, and worker orchestration.
- `physicsWorker.js` — numerical model and heatmap generation (`shear()` implements the ξ+ model
  used by this lab).
- `data/reference.json` — small, auditable survey-style reference-data bundle plotted alongside
  the simulation.
- `data/research-reference.json` — separate benchmark anchors used by the repository-quality
  validator, independent of the physics model.
- `research-overlay.js` — non-invasive quality/telemetry overlay.
- `scripts/validate.js`, `scripts/validate_repository.mjs` — no-dependency repository validation.

## Math Appendix

**Two-point shear correlation functions.** For a convergence power spectrum P_κ(ℓ), the shear
correlation functions are

```
ξ+(θ) = (1/2π) ∫ dℓ · ℓ · P_κ(ℓ) · J0(ℓθ)
ξ-(θ) = (1/2π) ∫ dℓ · ℓ · P_κ(ℓ) · J4(ℓθ)
```

where J0 and J4 are Bessel functions of the first kind. P_κ(ℓ) itself follows from a Limber
projection of the matter power spectrum P(k, z) weighted by the lensing efficiency kernel of the
source redshift distribution.

**Amplitude scaling (S8).** The overall amplitude of the shear correlation approximately factors
as

```
ξ+(θ) ∝ S8^2,   S8 = σ8 · sqrt(Ω_m / 0.3)
```

which is why cosmic shear surveys typically report constraints in the (σ8, Ω_m) plane compressed
along the S8 direction, rather than on σ8 alone.

**This lab's implementation.** `physicsWorker.js` uses a compact closed-form stand-in for the
Limber/Bessel integral above, calibrated to reproduce the qualitative amplitude and shape
scalings:

```
amp   = 2.5e-4 · (σ8 / 0.81)^2 · (Ω_m / 0.3)^1.4 · (z_source / 0.75)^0.8
ξ+(θ) = amp · (θ / 2)^-0.72 · exp(-θ / 220) · (1 - 0.12 · ia)
S8    = σ8 · sqrt(Ω_m / 0.3)
```

with θ in arcminutes. The (θ/2)^-0.72 term reproduces the roughly power-law decline of ξ+ at
small-to-intermediate angular scales, the exp(-θ/220) term suppresses the signal at large scales,
and the (1 - 0.12·ia) factor is a simplified linear model of intrinsic-alignment contamination
(non-lensing correlations in intrinsic galaxy shapes, which partially cancel or add to the true
lensing signal depending on sign and survey geometry).

## References

- Bartelmann, M. and Schneider, P., 2001. Weak gravitational lensing. Physics Reports, 340(4-5),
  pp.291-472.
- Troxel, M.A. et al., 2018. Dark Energy Survey Year 1 results: Cosmological constraints from
  cosmic shear. Physical Review D, 98(4), p.043528.
- Hildebrandt, H. et al., 2020. KiDS+VIKING-450 and DES-Y1 combined: Mitigating baryon feedback
  uncertainty with COSEBIs. Astronomy & Astrophysics, 638, A96.
- Asgari, M. et al., 2021. KiDS-1000 cosmology: Cosmic shear constraints and comparison between
  two point statistics. Astronomy & Astrophysics, 645, A104.
- Wilson, G. et al., 2017. Good enough practices in scientific computing. PLOS Computational
  Biology, 13(6), p.e1005510.

## The S8 Tension: Real Published Constraints

The main plot now uses **log-log axes** -- `xi_plus` spans roughly five decades in amplitude
over the plotted angular range, and a linear plot compresses nearly all of that structure
into a sliver near the y-axis. `data/reference.json` also carries a `published_S8_values`
block with real published `S8 = sigma8 * sqrt(Omega_m/0.3)` constraints for direct
comparison against the model's own `S8` telemetry metric:

- KiDS-1000 (Asgari et al. 2021): `S8 = 0.759 (+0.024/-0.021)`
- DES Y3 (Amon et al. 2022 / Secco et al. 2022): `S8 ~ 0.759 +/- 0.025`
- Planck 2018 CMB (Planck Collaboration 2020): `S8 = 0.834 +/- 0.016`

The mild ~2-3 sigma pull between low-redshift weak-lensing `S8` values and the Planck CMB
value is known in the literature as the **S8 tension** -- drag the `sigma8`/`omegaM` sliders
until the live `S8` metric lands near 0.76 versus 0.83 to build intuition for how large a
shift that actually is.

## Research Quality Upgrade

See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for the validation layer, reference anchors,
equations, and research boundaries added to this repository.
