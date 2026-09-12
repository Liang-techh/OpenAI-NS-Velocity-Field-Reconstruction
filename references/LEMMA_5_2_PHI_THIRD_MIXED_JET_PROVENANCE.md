# Lemma 5.2 repaired `phi_n` third-mixed-jet provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: 50-52.

Implemented identities:

- Lemma 5.2 / Eq. (5.14) repairs the angular coefficient by compact bumps,
  `E_n = E_tilde_n + sum_j beta_{n,j}(eta) b^E_j(R)`.
- The PositiveAxis variable is recovered by the exact relation `phi_n = C E_n / R`, with `R=sqrt(2X)` and the smooth axis extension inherited from the inner solution.
- Therefore `Delta phi_n = C sum_j beta_{n,j}(eta) b^E_j(R)/R` on the compact repair supports.
- The preceding-diffusion eta derivative used by the strict-lower source requires the scalar third-mixed fields `phi_XXeta`, `phi_Xetaeta`, and `phi_etaetaeta`.

Writing `g_j(X)=b^E_j(sqrt(2X))/sqrt(2X)`, the executable bridge adds

- `C sum_j beta'_j g_j,XX` to `phi_XXeta`,
- `C sum_j beta''_j g_j,X` to `phi_Xetaeta`, and
- `C sum_j beta'''_j g_j` to `phi_etaetaeta`.

The landed `background_moment_repair_phi_jets.py` formulas for `g_j`, `g_j,X`, and `g_j,XX` are reused exactly. Every E-repair bump is supported on a strictly positive radial interval, so `1/R` is evaluated only for `R>0`; off support, including the axis, the compact correction is identically flat.

## Executable mapping

`src/openai_ns_reconstruction/background_moment_repair_phi_third_mixed_jets.py` provides `Lemma52RepairedPhiThirdMixedJetAdapter`.

The adapter takes:

1. the landed `Lemma52RepairedProfileAdapter`, which owns the actual compact Lemma-5.2 repair geometry and ordinary eta derivatives of the repair coefficients;
2. the manuscript normalization `C`;
3. a caller-supplied unrepaired/inner `ProfileThirdMixedJet` for `phi_n`.

On active repair support it requests repair-coefficient derivatives through order three and returns one coherent `ProfileThirdMixedJet`. Off support it returns the supplied base jet exactly and does not request unnecessary high-order moment/patch data. No production finite differences are used.

`tests/test_background_moment_repair_phi_third_mixed_jets.py` checks that the lower six components exactly equal the pre-existing repaired-phi second-jet bridge. It then independently centered-differences that older second-jet path in eta and compares the three new mixed fields. The test also covers exact axis/off-support preservation and fail-closed behavior when third eta repair data are unavailable on active support.

## Truth boundary

This increment is **formal-structure / solver infrastructure**, not a paper-exact coefficient construction.

The unrepaired/inner `phi_n` third-mixed jet, normalization `C`, unrepaired moments, patch factor, and their eta derivatives remain upstream inputs. In particular this increment does not itself make `partial_eta actualLowerSource` hierarchy-owned and does not materialize the genuine `k=1` Eq. (5.7) Picard step. Issue #1 still owns the real leading profile data.

The common analytic strip, coefficient bounds (5.17), converged all-order PositiveAxis solve, recursive cutoff/local-finiteness argument, Proposition 5.3 all-jets residual decay, and paper-exact velocity remain open. No paper-exact gate is changed and `paper_exact_velocity_available` remains false.
