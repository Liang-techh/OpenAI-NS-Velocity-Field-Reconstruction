# Lemma 5.2 compact repaired phi fifth-mixed jet provenance

## Scope

This increment lifts the already-landed compact Lemma 5.2 swirl repair from a
fourth-mixed `phi_n` jet to the next mixed derivative tier needed on the route
to a hierarchy-owned third-eta strict-lower source.

On active E-repair supports,

`Delta phi_n = C sum_j beta_j(eta) b_j(R) / R`, with `R=sqrt(2X)`.

Because the radial bump factors are eta-independent, the three newly required
total-order-five entries are obtained analytically from the same repair
coefficient jet:

- `phi_XXetaetaeta` from `beta''' * d_X^2(b/R)`;
- `phi_Xetaetaetaeta` from `beta'''' * d_X(b/R)`;
- `phi_etaetaetaetaeta` from `beta''''' * (b/R)`.

All lower fields are delegated to the previously landed repaired fourth-mixed
bridge using the exact fourth-order projection of the same base fifth-mixed
provider.  Off compact repair support the base fifth-mixed provider is returned
exactly and no unnecessary fifth eta row is requested.

## Paper / implementation map

- Lemma 5.2 and Appendix A: compact five-bump moment repair.
- Eqs. (5.3)-(5.6): strict-lower recursive source whose higher eta jets require
  corresponding higher mixed profile jets.
- Eq. (5.7): downstream PositiveAxis Picard solve.  The landed exact
  second-eta Picard rule exposes a third-eta previous-iterate dependency, which
  in turn requires a hierarchy-owned third-eta forcing/source path.

The implementation reuses the exact compact `b/R` X-second-jet helper,
`Lemma52RepairedProfileAdapter.repair_coefficients(max_order=5)`, and the
existing fourth-mixed phi repair.  It does not introduce another moment solver,
generic cutoff, sampled radial differentiation, or production finite
difference.

## Verification

Regression uses the real `Lemma52MomentRepair.from_intervals(...)` five-bump
repair with analytic degree-five moment and patch-factor jets.

1. The complete repaired fifth-mixed jet projects exactly to an independently
   instantiated landed repaired fourth-mixed bridge.
2. The three new total-order-five fields are compared against centered eta
   differences of the old analytic fourth-mixed path on active compact support.
   Finite difference is test oracle only.
3. Outside compact support, incomplete higher moment/patch jets are not touched
   and the base fifth-mixed jet is returned exactly; the same incomplete data
   fail closed on active support.
4. Invalid normalization and non-jet base providers are rejected.

## Truth boundary

Status remains **Stage 2 / `formal-structure`** with
`full_reconstruction=false` and `paper_exact_velocity_available=false`.

The unrepaired/base fifth-mixed phi jet and the moment/patch eta jets remain
explicit upstream contracts.  Issue #1 still must supply the genuine strong
leading profile before this derivative layer can participate in a paper-backed
recursive solve.  Caller-supplied jets, polynomial fixtures, generic cutoffs,
and finite-order numerical fits are not paper-exact.

This increment does not make the fifth-mixed phi provider hierarchy-owned and
does not yet implement a third-eta angular preceding-diffusion operator or a
hierarchy-owned `partial_eta^3 actualLowerSource`.  Consequently it does not
claim hierarchy-owned `partial_eta^3 f_n`, `partial_eta^2 W_n^(1)`, Picard
convergence, final positive-order coefficient materialization, a paper recursive
cutoff-scale schedule, Proposition 5.3 all-jets residual decay, or a paper-exact
reconstructed velocity.
