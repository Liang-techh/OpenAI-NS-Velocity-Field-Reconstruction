# Section 7 large-band phase/frame bridge provenance

Status: **formal-structure**. This bridge does not make the velocity paper-exact.

## Paper / pinned-source mapping

`src/openai_ns_reconstruction/phase_large_band_frame.py` connects two already-landed theorem-side adapters at one integer dyadic band.

- Section 6 / Section 7.1: `S_*=ell^2`, `epsilon=2^(-ell h)` and the rounded carrier scale feed the phase-normal estimates leading to Eqs. (7.9)-(7.11).
- `PhaseEstimates.phaseError_le_four_div`: the existing `LargeBandPhaseScaleCertificate` proves the scalar large-band hypotheses without constructing binary64 `Q=2^-ell` and yields the local consequence `4*phaseConstant(M)/S_*`.
- Pinned `BasePhaseGeometry.frame_errors_of_normal_close` and `damping_error_of_normal_close`: the existing `FrameDampingEnvelope.from_large_band` uses the family-level specialization `A=3M`, `delta=BasePhaseGeometry.phaseConstant(M)/S_*`, `eta=16 M^2/S_*`, `b=normalLower(M,u)`, and checks `b<=B`, `delta<=B/2`, `B<=M`, `0<=nu<=4` before exposing the frame and damping envelopes.

The implementation deliberately keeps the local `PhaseEstimates.phaseConstant` bound and the enlarged family-level `BasePhaseGeometry.phaseConstant` bound separate. They are different constants in the pinned source and are not identified by this bridge.

## Independent regression

`tests/test_phase_large_band_frame.py` rebuilds the family constants directly from their closed formulas rather than calling the production constant helpers. For the explicit formal-structure fixture `h=1/200`, `M=2`, `u=1`, `B=1`, it independently checks that the exact phase-scale gate is already valid while the family normal-smallness inequality fails at `ell=94118` and first passes at `ell=94119`. It separately reconstructs the normal, coordinate-frame, and damping constants and compares the resulting `1/S_*` envelopes.

The fixture is theorem arithmetic only; it is not paper field data and is never used as a surrogate wave.

## Remaining boundary

Passing this bridge certifies only scalar compatibility of the large-band phase gate and the frame/damping implication constants. It does **not** certify the Proposition 5.5 base field, the `LocalBaseBounds` C1/C2 hypotheses on any active slow box, the vector normal/slot-derivative hypotheses, the actual viscosity/base comparison coming from the manuscript field, support separation, a stress decomposition, an amplitude solution, a curl wave, or uniform Eqs. (7.9)-(7.11) for the paper construction. `paper_exact_velocity_available` therefore remains `false`.
