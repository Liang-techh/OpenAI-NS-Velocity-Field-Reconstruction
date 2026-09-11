# Lemma 5.2 repaired `phi_n` second-jet provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: 50-52.

Implemented identities:

- Lemma 5.2 / Eq. (5.14) repairs the angular coefficient by compact bumps,
  `E_n = E_tilde_n + sum_j beta_{n,j}(eta) b^E_j(R)`.
- Immediately after the repair the PositiveAxis variable is recovered by the exact relation
  `phi_n = C E_n / R`, with `R=sqrt(2X)` and the smooth axis extension inherited from the inner solution.
- Hence the compact correction is
  `Delta phi_n = C sum_j beta_{n,j}(eta) b^E_j(R)/R`.
- Because every repair bump is supported on a strictly positive radial interval, the `1/R` factor is differentiated only where `R>0`. Off the repair supports the compact correction and all of its derivatives vanish exactly.

For `g(X)=b(sqrt(2X))/sqrt(2X)`, the executable bridge uses

- `g_X = b_X/R - b/R^3`,
- `g_XX = b_XX/R - 2 b_X/R^3 + 3 b/R^5`,

and combines these radial derivatives with the already-landed ordinary eta derivatives of the repair coefficients through order two.

## Executable mapping

`src/openai_ns_reconstruction/background_moment_repair_phi_jets.py` provides `Lemma52RepairedPhiSecondJetAdapter`.

The adapter takes:

1. the landed `Lemma52RepairedProfileAdapter`, which owns the actual compact Lemma-5.2 repair geometry and eta-dependent repair coefficients;
2. the manuscript normalization `C`;
3. a second-jet provider for the unrepaired/inner `phi_n` profile.

It returns a `ProfileSecondJet` for the repaired `phi_n`, suitable for the landed strict-lower-history PositiveAxis source interface. No production finite differences are used. At `X=0` and everywhere outside the compact E-repair supports, the supplied inner/base `phi_n` jet is returned exactly and no `C/R` division is evaluated.

`tests/test_background_moment_repair_phi_jets.py` independently differentiates the function-level identity

`phi_repaired = phi_base + C (E_repaired-E_base)/sqrt(2X)`

with centered finite differences on an active compact support and compares all six second-jet components. It also verifies exact off-support/axis preservation and fail-closed behavior when the active repair lacks the required eta-jet rows.

## Truth boundary

This increment is **formal-structure / solver infrastructure**, not a paper-exact coefficient construction.

The unrepaired inner `phi_n` second jet is still an upstream input. The normalization `C`, unrepaired moments, patch factor, and their eta derivatives are also supplied rather than derived from a completed recursive hierarchy. This module therefore does not establish the common analytic strip, coefficient bounds (5.17), the full converged positive-order solve, the infinite recursive cutoff/local-finiteness argument, or Proposition 5.3 all-jets-flat residual decay.

No paper-exact gate is changed and `paper_exact_velocity_available` remains false.
