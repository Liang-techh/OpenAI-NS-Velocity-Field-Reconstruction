# Section 7 phase-estimate provenance

Source pin: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Reviewed formalization module: `NavierStokes/PhaseEstimates.lean`.

This repository's `src/openai_ns_reconstruction/phase_estimates.py` transcribes only the scalar scale-control layer used by the formalized rounded phase estimates:

- `phaseError(S, epsilon, k) = 1/S + S*epsilon^2 + S/k + epsilon*S`;
- `phaseConstant(M) = 8*M^3 + 2*M^4`;
- theorem `phaseError_le_four_div`, whose hypotheses are `S>0`, `S^2*epsilon^2<=1`, `S^2/k<=1`, and `epsilon*S^2<=1`, and whose conclusion is `phaseError<=4/S`;
- the downstream `rounded_normal_estimates` bound uses `phaseConstant(M)*phaseError` for both the normal error and normal-velocity bound once its additional representative-frequency and local-base derivative hypotheses are supplied.

The executable `PhaseScaleCertificate` intentionally checks only these scale hypotheses plus the downstream domain assumptions `M>=1`, `S>=1`, `k>=1`, `epsilon>=0`. It does **not** claim that the paper's slow boxes have been selected, that the local C2/base-field hypotheses hold, that Eqs. (7.9)-(7.11) are uniform over all active boxes, or that any oscillatory wave is paper-exact. Those remain Issue #3 blockers.
