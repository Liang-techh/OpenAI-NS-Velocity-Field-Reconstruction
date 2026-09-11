# Section 7 frame/damping implication provenance

Status: **formal-structure**. This file does not certify a paper-exact base field or oscillatory wave.

Pinned formal source:

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Module: `NavierStokes/BasePhaseGeometry.lean`

Implemented in `src/openai_ns_reconstruction/phase_frame_bounds.py`:

- `damping_denominator` ↔ `BasePhaseGeometry.dampingDenominator`.
- `normal_lower` ↔ `BasePhaseGeometry.normalLower`.
- `normal_constant` ↔ `BasePhaseGeometry.normalConstant M = 8 * PhaseEstimates.phaseConstant (2*M)`.
- `frequency_bound` ↔ the pinned `BasePhaseGeometry.frequencyBound` enlargement that uniformly dominates the frozen phase/frequency constants.
- `base_phase_constant` ↔ the later family-level definition `BasePhaseGeometry.phaseConstant M = normalConstant (frequencyBound M)`. This distinction is essential: the family-level constant is deliberately larger than `normalConstant M`.
- `frame_coordinate_error_bound` ↔ the conclusion of `frame_errors_of_normal_close`: with `G=M+2+2A` and `E=eta+8(1+A)delta/b`, each of the three changing-frame coefficient discrepancies is bounded by `16 G^2 (1+G) E`.
- `damping_error_bound` ↔ the conclusion of `damping_error_of_normal_close`: `|nu*||n||^2 - nu*B^2(1+s^2)| <= 4 M (2A+5) delta` for `0<=nu<=4` under the theorem's normal-closeness hypotheses.
- `coordinate_constant` and `damping_constant` ↔ the named family constants used by `FamilyData.coordinate_errors`, `FamilyData.damping_error`, and `FamilyData.coefficientControl` after the specialization `A=3M`, `delta=phaseConstant(M)/S`, `eta=16 M^2/S`, `b=normalLower(M,u)`.
- `FrameDampingEnvelope.from_large_band` checks the scalar hypotheses `normalLower<=B<=M`, `delta<=B/2`, and `0<=viscosity<=4` before exposing the specialized error envelopes.

Independent checks in `tests/test_phase_frame_bounds.py` separately recompute `normalConstant`, `frequencyBound`, and the enlarged family `phaseConstant`; recompute the frame/damping theorem arithmetic directly; verify the large-band specialization reduces exactly to `coordinateConstant/S` and `dampingConstant/S`; and exercise fail-closed paths for every scalar smallness/boundedness hypothesis.

## What this does not prove

The official theorems also require vector/analytic facts that cannot be inferred from scalar samples: `||K||=1`, normal closeness, slot-normal derivative closeness, representative shear orthogonality, actual/reference `F` and shear comparisons, and the differentiability/convexity assumptions feeding `LocalBaseBounds`. Those remain upstream obligations and are not encoded as guessed booleans or fitted data here.

Consequently this module does **not** instantiate Proposition 5.5 base fields, does not complete the uniform Eqs. (7.9)–(7.11) proof on every active slow box, does not solve the amplitude ODE, and does not construct a stress-cone/curl wave. It must remain `formal-structure`; `paper_exact_velocity_available` stays false until the actual background and all downstream construction layers are materialized.
