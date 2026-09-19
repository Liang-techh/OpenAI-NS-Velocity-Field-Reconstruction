# Validation protocol and current status

Three questions must stay separate: whether files/code are intact, whether finite-sample numerical gates pass, and whether a continuous mathematical assertion has been proved. A green workflow answers only the tasks actually executed inside it.

## Original full protocol

The independently implemented Cartesian operator is `_runtime/validate.py`, imported with unchanged formulas. For each seed it uses 4,096 uniform points in `[-2,2]^3` and times `[0.25,0.3125,0.4375,0.5625,0.6875,0.75]`.

Spatial refinement uses steps `0.02,0.01,0.005` at fixed time step `0.0025`. Time derivatives use a centered fourth-order stencil where possible and one-sided fourth-order stencils at the window boundaries. A separate 256-point subset varies time steps `0.02,0.01,0.005` at fixed spatial step. Energy quadrature independently uses orders 24, 48 and 96.

The residual norms at each time are

```text
sampled maximum = max_i ||R_i||_2
spatial volume L2 estimate = sqrt(64 * mean_i ||R_i||_2^2)
```

Both original momentum gates require values below `0.001`. Other gates check divergence, initial energy, energy range/convergence, sampled boundary zeros, signs at an original core probe, scaled-probe drift and the derived velocity amplitude. These other gates do not replace the failed momentum gates.

## Replays, not new holdouts or new optimization

The original ST054 study froze its fields before seeds 9175491 and 9175492. This release reuses those seeds to test reproducibility. It does not pretend they are newly independent of all prior research. The source summary is preserved unchanged under `evidence/upstream_ST054_results.json`.

Run `python tools/replay_release.py --out outputs/full_release_replay` to repeat four complete validations. Each scientific result is a rejection. The replay script returns success only when the known rejected gate booleans and reference norms are reproduced within the stated `1e-5` cross-platform comparison tolerance. That tolerance is a comparison tolerance, not a new acceptance threshold.

Individual validation uses `ns-field validate ...` and returns `1` when scientific gates fail. Existing output files are never overwritten. Reports identify the exact release candidate JSON hash. The original research JSON hash is a separate provenance field; regenerated metadata and tiny upstream reconstruction rounding are not misrepresented as byte identity.

## What the plots do not measure

A displayed grid maximum can miss peaks. The original denser-grid study still found maxima near `0.038` for these candidates. Neither a sampled maximum nor a finite-difference convergence trend provides a continuous-space supremum certificate. Local geometry statistics depend on observation window, weighting, seeding, arc length and resolution. Direction fractions apply only to stated finite probes.

The source study also identifies an axisymmetric-pressure limitation: fixed velocity and force make the angular residual independent of pressure. Its nonzero measured angular component is one reason pressure-only optimization cannot finish the target. Future work must change the actual velocity dynamics, not a color scale or threshold.
