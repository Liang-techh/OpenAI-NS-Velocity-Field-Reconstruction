# Repaired-history phi third-mixed provenance

## Scope

This increment connects the already-landed analytic Lemma 5.2 compact
`phi_n = C E_n / R` third-mixed jet to a strict-lower Section 5 coefficient
history without changing the older compatibility hierarchy.  It is solver
infrastructure for the eta derivative of the Eq. (5.3)-(5.6) lower-order source.

The implementation is
`src/openai_ns_reconstruction/background_repaired_history_phi_third_mixed.py`.
It introduces a strong coefficient-source subtype that owns one
`Lemma52RepairedPhiThirdMixedJetAdapter` together with the already-landed
fourth-mixed repaired `U_n` adapter.  The inherited phi second jet is not a
second caller-maintained record: it is the exact `ProfileThirdMixedJet.second()`
projection of that same repaired provider.  The inherited U third-mixed path is
likewise the projection of the same fourth-mixed U provider, and beta remains
derived only by the analytic Eq. (5.2) adapters in the base hierarchy.

## Paper connection

The relevant paper layer is Section 5, especially Eqs. (5.3)-(5.7), Lemma 5.2,
and the compact moment correction described around Eq. (5.14) and Appendix A.
The immediate downstream use is the analytic eta derivative of
`precedingDiffusion = Z_(b-D)(Z_b F_(n-1))`, which requires
`phi_XXeta`, `phi_Xetaeta`, and `phi_etaetaeta`.

No new formula for the Lemma 5.2 correction is invented here.  The strong
history source delegates those derivatives to the previously landed
`Lemma52RepairedPhiThirdMixedJetAdapter` and only establishes coherent ownership
inside the recursive coefficient history.

## Fail-closed boundary

The Lemma 5.2 factory is positive-order only.  A strong order-zero source can be
constructed explicitly, but this module does not fabricate it: the actual
leading third-mixed profile remains an Issue #1 input.  Therefore an attempted
strong phi query on an older source fails with a missing-data error.  A source
whose third-mixed jet does not project exactly to the hierarchy-owned second jet
also fails closed.

The following remain upstream or incomplete:

- the genuine Issue #1 leading third-mixed profile;
- unrepaired positive-order base phi/U mixed jets and the moment/patch eta-jets;
- a complete hierarchy-owned eta derivative of `actualLowerSource` and `f_n`;
- the genuine hierarchy-owned `k=1` Eq. (5.7) Picard application;
- recursive coefficient materialization, recursive cutoff-scale completion, and
  Proposition 5.3 all-jets truncation/residual decay;
- paper-exact background or final velocity.

Accordingly Stage 2 remains `formal-structure` and
`paper_exact_velocity_available=false`.

## Validation

`tests/test_background_repaired_history_phi_third_mixed.py` uses the actual
five-bump `Lemma52MomentRepair` geometry and polynomial eta moment jets.  It
cross-checks the hierarchy-owned repaired phi third-mixed jet against an
independently instantiated `Lemma52RepairedPhiThirdMixedJetAdapter`, checks that
the ordinary phi second jet is exactly its projection, and verifies that the
same strong coefficient source still owns the fourth-mixed repaired U layer.
It also verifies fail-closed behavior for an older fourth-mixed-U source without
third-mixed phi ownership and for a deliberately incoherent third-to-second phi
projection.

These are floating-point regression checks of the implemented formulas and
ownership invariants.  They are not a proof of the paper's all-order estimates,
uniform moment invertibility, recursive cutoff schedule, or residual decay.
