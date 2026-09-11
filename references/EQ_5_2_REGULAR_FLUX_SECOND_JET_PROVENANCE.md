# Eq. (5.2) regular-flux second-jet provenance

Status: **formal-structure only**.

This increment does not materialize a recursively solved Section-5 coefficient
and does not change `paper_exact_velocity_available=false`.

## Pinned source and dependency chain

Official Lean source commit:
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Paper locations: Eq. (5.2), used downstream by Eqs. (5.3)-(5.7) and the
strict-lower PositiveAxis recursion.

The repository already had:

- `background.py::coefficient_radial_flux`, implementing the Eq. (5.2)
  reconstruction of `V_n` from a coefficient profile;
- `background_lower_history_source.py::ProfileSecondJet`, whose `beta_jets`
  entries represent the regular quotient `beta_j=V_j/X` needed by the exact
  lower-history source and Eq. (5.6) evaluator; and
- `background_moment_repair_second_jets.py`, which supplies repaired `U_n` and
  `E_n` second jets but explicitly stopped short of constructing `beta_n`.

## What this increment adds

`background_regular_flux_second_jets.py` differentiates Eq. (5.2) analytically
and returns the full second `(X,eta)` jet of `beta_n=V_n/X` for
`lambda_n=2 n h`.

Writing

`A[U](X,eta)=integral_0^1 U(sX,eta) ds`,

`D=1/2-h`, `d=1-eta^2`, and `L=1-2h eta^2`, the implemented value is

`beta_n = [2 eta U_n - 2 (D+lambda_n) eta A[U_n] - d A[partial_eta U_n]] / L`.

The radial-average derivatives are evaluated using the exact chain-rule identity

`partial_X^a partial_eta^b A[U](X,eta)
 = integral_0^1 s^a partial_X^a partial_eta^b U(sX,eta) ds`.

Consequently no production path samples `V_n` and divides by `X`; the axis
`X=0` is handled directly by the regular formula.

A full second jet of `beta_n` requires three total-order-three derivatives of
`U_n`: `U_XXeta`, `U_Xetaeta`, and `U_etaetaeta`.  The new
`AxialThirdMixedJet` type makes this requirement explicit rather than silently
estimating those derivatives from finite differences or lower-order data.

The radial integrals use the repository's cached Gauss-Legendre rules.  This is
an executable numerical realization of the exact integral identities, not a
rigorous quadrature certificate.

## Independent regression checks

`tests/test_background_regular_flux_second_jets.py` uses a separable polynomial
`U(X,eta)=P(X)Q(eta)` for which `A[U]` is integrated independently in closed
form.  The test oracle evaluates the resulting closed-form Eq. (5.2) value and
recovers its first, pure-second, and mixed derivatives with centered finite
stencils.  Those numerical derivatives are compared against the production
analytic second jet.

Additional tests verify that:

1. `X=0` is finite and agrees with the closed-form regular quotient without any
   `V/X` division;
2. missing total-order-three jet data fail closed rather than being inferred;
3. invalid coefficient order and domain inputs fail closed.

Finite differences are test-only oracles; production differentiation is
analytic.

## Truth boundary / remaining blocker

The new adapter still consumes an `AxialThirdMixedJet` provider.  It therefore
does **not** prove that the supplied `U_n` derivatives come from the converged
Eq. (5.7) positive-order coefficient or from the Lemma-5.2 repaired hierarchy.
In particular, the currently landed repaired-profile bridge only supplies
second jets; hierarchy-derived third mixed eta derivatives and the corresponding
third-order repair-coefficient jets must still be materialized before this
adapter can automatically populate the strict-lower `beta_jets` history.

This increment also does not prove a common analytic strip, construct the full
positive-order Picard fixed point, certify the infinite recursive cutoff
schedule/local finiteness, or establish Proposition 5.3 all-jets-flat residual
decay.  Stage 2 therefore remains `formal-structure`.
