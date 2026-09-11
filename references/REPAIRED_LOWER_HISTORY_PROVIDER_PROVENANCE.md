# Section 5 repaired lower-history provider provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: 45-52.

This increment does not add a new manuscript identity. It connects previously landed paper/Lean-exact interfaces that were still manually separated:

- Eqs. (5.3)-(5.7): the PositiveAxis recursion consumes strict lower history in the variables `phi_j`, `U_j`, and the regular quotient `beta_j=V_j/X`;
- Eq. (5.2): `beta_j` is determined by the axial coefficient `U_j`, so it must not be supplied independently when a coherent coefficient is already owned;
- Lemma 5.2 / Eq. (5.14): positive-order `U_j` and `E_j` are compactly repaired, with `phi_j=C E_j/R` on the positive-radius repair support;
- the already-landed `actualLowerSource` / `lowerHistoryData` formulas use only orders strictly below the current recursive order.

## Executable mapping

`src/openai_ns_reconstruction/background_repaired_history.py` adds two fail-closed containers.

`Section5CoefficientJetSource` owns one coefficient's PositiveAxis `phi` second-jet provider and one coherent axial `U` third-mixed-jet provider. Its `from_lemma52_repair(...)` constructor instantiates the landed `Lemma52RepairedPhiSecondJetAdapter` and `Lemma52RepairedThirdMixedJetAdapter` for a positive order using the actual compact Lemma-5.2 repair geometry. It explicitly rejects order zero because the leading profile belongs to Issue #1 and is not a Lemma-5.2 positive-order repair.

`Section5LowerHistoryJetHierarchy` owns a contiguous sequence of coefficient sources `0,...,N`. It:

1. projects each owned coherent axial third-mixed jet to the second jet required by `actualLowerSource`;
2. derives `beta_j=V_j/X` only through the landed analytic Eq. (5.2) adapter, so callers cannot pair an owned `U_j` with an unrelated beta jet;
3. rejects missing/gapped orders and rejects Lemma-5.2 sources repaired with a different normalization `C`;
4. exposes `positive_axis_providers(n)` and `positive_axis_fields(n)`, which feed the hierarchy directly into the landed strict-lower-history / Eq. (5.7) bridge and query only `0,...,n-1`;
5. exposes `first_picard_term(n,xi,eta)` as the already-landed genuine `G f_n` step driven only by the owned strict lower history.

`tests/test_background_repaired_history.py` uses the repository's explicit analytic toy profile only as a test oracle. The positive-order entry is nevertheless repaired by the actual `Lemma52MomentRepair` / `Lemma52RepairedProfileAdapter`, not by a pointwise toy correction. The tests independently compare the hierarchy-owned repaired phi/U jets with the landed analytic repair adapters, compare the hierarchy-derived Eq. (5.2) beta value with the function-level repaired profile's `radial_flux_factor`, verify that order `n=2` can build actual lower-source / Eq. (5.7) fields while no order-2 coefficient exists, and verify fail-closed behavior for order gaps and normalization mismatch.

## Truth boundary

This increment is **formal-structure / solver infrastructure**, not a paper-exact coefficient solve.

The order-zero source is still an explicit upstream provider because Issue #1 has not yet supplied the final materialized Theorem 4.6 profile to this hierarchy. Positive-order base jets, unrepaired moments, patch-factor jets, and their eta derivatives also remain inputs to the already-landed Lemma 5.2 repair adapters. The new hierarchy therefore removes a manual *wiring* gap but does not establish existence or convergence of the full recursive coefficient family.

In particular, this increment does **not** prove a common analytic strip, coefficient bounds (5.17), an all-order recursive cutoff/local-finiteness schedule, Proposition 5.3 all-jets residual decay, or a completed cutoff-summed paper-exact background. No completion gate is changed and `paper_exact_velocity_available` remains false.
