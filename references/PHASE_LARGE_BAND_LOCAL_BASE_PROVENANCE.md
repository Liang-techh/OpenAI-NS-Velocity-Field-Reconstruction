# Large-band LocalBase admission provenance

Status: **formal-structure**. This increment does not materialize a paper-exact background and does not upgrade Stage 3 truth status.

## Why this layer exists

The executable Section 6 bridge in `slow_labels.py` intentionally relies on the binary64 `DyadicChart`, which is capped at `ell<=1000`. The uniform Section 7 phase/frame estimates are eventual large-band statements, and the exact gate in `phase_large_band_scale.py` can certify much larger integer bands without evaluating `Q=2^-ell`.

At those bands it is also undesirable to encode the `LocalBaseBounds` hypotheses as raw binary64 errors such as `M epsilon^2`: eventually `Q`, `epsilon`, or `epsilon^2` can underflow even though the analytic theorem remains meaningful.

`src/openai_ns_reconstruction/phase_large_band_local_base.py` therefore adds a theorem-facing, underflow-safe admission layer. It does **not** replace the executable slow-label geometry and it does not infer any field bound from samples.

## Paper / Lean mapping

Pinned upstream source: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, primarily `NavierStokes/PhaseEstimates.lean`, together with the Section 6 label `gamma=(ell,a,sigma)` and mesh scale `S_*^-3`.

The normalized witness fields encode the quantitative content of `PhaseEstimates.LocalBaseBounds`:

- second derivative bounds for `F0,G0`, first derivative bounds for `F0`, and directional bounds for the actual `F,G` fields are divided by `M` and must be `<=1`;
- actual/reference C0/C1 errors are divided by `M epsilon^2` and must be `<=1`;
- the certified slow-domain diameter is stored as `d S_*^3` and must be `<=1`, exactly representing `d<=S_*^-3`;
- convexity, differentiability, representative membership, and the uniform-supremum theorem obligation must be separately certified true.

The admission is tied to the exact `LargeBandPhaseScaleCertificate` by the same integer `ell` and exactly the same `M`. It exposes `log2(epsilon^2)=-2 ell h` as a `Fraction` and `S_*^-3=ell^-6` exactly, so no binary64 `Q` or `epsilon` is required.

## Evidence boundary

Only evidence tagged `analytic-theorem` or `formal-theorem` is admitted, and a nonempty provenance description is mandatory. Tags such as sampled/fitted/numeric-scan are rejected. This is a fail-closed interface convention, not a machine checker for the cited theorem itself.

`AsymptoticSlowLabel` is only a large-band label address. It does not claim that a `DyadicChart`, partition cutoff, slow representative, or support has been numerically materialized. Likewise, admission of normalized LocalBase hypotheses does not prove that the underlying Proposition 5.5 fields are the manuscript fields.

## Independent regression

`tests/test_phase_large_band_local_base.py` uses `ell=100000`, `h=1/200` and independently verifies:

- binary64 `2^-ell` has underflowed while the admission remains usable;
- `S_*=ell^2`, `S_*^-3=ell^-6`, and `log2(epsilon^2)=-1000` are recovered exactly;
- all normalized LocalBase thresholds fail closed above one;
- sampled/fitted evidence, missing provenance, missing qualitative theorem facts, mismatched band identity, and mismatched `M` are rejected;
- sign-paired labels retain the same sign-free slow-box key with no `ell<=1000` cap.

The test values are formal fixtures only and are not samples from the paper's flow.

## Remaining boundary

Still missing for paper-exact Eqs. (7.9)--(7.11):

- the actual Proposition 5.5 paper background and a machine-linked theorem proving these LocalBase ratios for every active enlarged slow box;
- the actual smooth Section 6 squared partition, representatives, common rectangles and support separation at asymptotic bands;
- field-level vector/normal hypotheses needed by the frame and damping theorems;
- stress-cone amplitudes, curl waves, mean corrections, and the Section 9 correction sequence built from those actual fields.

Therefore `actual_base_fields_verified=false`, `uniform_eq_7_9_to_7_11_verified=false`, and `paper_exact_velocity_available=false` remain mandatory.
