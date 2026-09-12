# Lemma 5.2 compact repaired phi fourth-mixed jet provenance

## Scope

This increment lifts the already-landed compact Lemma 5.2 swirl repair from a
third-mixed `phi_n` jet to the `ProfileFourthMixedJet` required by the analytic
second-eta preceding-diffusion operator.

On active E-repair supports,

`Delta phi_n = C sum_j beta_j(eta) b_j(R) / R`, with `R=sqrt(2X)`.

Because the radial bump factors are eta-independent, the three newly required
entries are obtained analytically from the same repair coefficient jet:

- `phi_XXetaeta` from `beta'' * d_X^2(b/R)`;
- `phi_Xetaetaeta` from `beta''' * d_X(b/R)`;
- `phi_etaetaetaeta` from `beta'''' * (b/R)`.

All lower fields are evaluated from the same formula and the resulting
fourth-mixed jet projects exactly to the previously landed repaired
third-mixed bridge.  Off compact repair support the base fourth-mixed provider
is returned exactly and no unnecessary fourth eta row is requested.

## Paper / implementation map

- Lemma 5.2 and Appendix A: compact five-bump moment repair.
- Eqs. (5.3)-(5.6): strict-lower recursive source.
- The angular preceding-diffusion term uses the same analytic
  `Z_(b-D)(Z_b F_(n-1))` implementation already landed for its second eta jet.
- Eq. (5.7): downstream PositiveAxis solve that will consume the eventual
  hierarchy-owned second forcing derivative.

The implementation reuses the existing exact compact `b/R` X-second-jet helper
and `Lemmas52RepairedProfileAdapter.repair_coefficients`; it does not introduce
a second moment solver, generic cutoff, or sampled radial differentiation.

## Verification

Regression uses the real `Lemma52MomentRepair.from_intervals(...)` five-bump
repair with analytic degree-four moment and patch-factor jets.

1. The complete repaired fourth-mixed jet projects exactly to an independently
   instantiated landed repaired third-mixed bridge.
2. The three new total-order-four fields are compared against centered eta
   differences of the old analytic third-mixed value path on active compact
   support.  Finite difference is test oracle only.
3. The repaired fourth-mixed phi jet is passed through the landed analytic
   second-eta preceding-diffusion operator, and that result is cross-checked
   against centered eta differentiation of the landed first-eta operator.
4. Outside compact support, incomplete higher moment/patch jets are not touched
   and the base fourth-mixed jet is returned exactly; the same incomplete data
   fail closed on active support.
5. Invalid normalization and non-jet base providers are rejected.

## Truth boundary

Status remains **Stage 2 / `formal-structure`** and
`paper_exact_velocity_available=false`.

The unrepaired/base fourth-mixed phi jet and the moment/patch eta jets are still
explicit upstream contracts.  In particular, Issue #1 still must supply the
genuine order-zero strong profile before this derivative layer can participate
in a paper-backed recursive solve.  Caller-supplied jets, polynomial fixtures,
generic cutoffs, and finite-order numerical fits are not paper-exact.

This increment does not yet make the fourth-mixed phi provider hierarchy-owned;
it therefore does not claim a hierarchy-owned `partial_eta^2 actualLowerSource`
or forcing, `partial_eta W_n^(1)`, Picard convergence, final positive-order
coefficient materialization, a recursive cutoff-scale schedule, Proposition 5.3
all-jets residual decay, or a paper-exact reconstructed velocity.
