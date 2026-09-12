# Section 9 weighted correction-jet product-rule provenance

Status: **formal-structure / analytic derivative adapter only**. This increment
does not materialize the manuscript's actual Section 7/8 correction fields,
does not certify the manuscript cutoff jet, and does not make paper-exact
velocity or forcing available.

## Scope

`section9_weighted_correction_jet.py` narrows one gap in the existing Eq. (9.21)
finite-prefix jet path. Instead of allowing a provider to supply derivatives of
the already-weighted quantity

`chi(a_j q) (A_j, B_j, p_j)`

directly, the new adapter takes two complete physical `(t,x,y,z)` derivative
jets through one total order:

1. the unweighted correction `(A_j,B_j,p_j)`; and
2. the scalar weight `chi(a_j q)`.

For every multi-index `alpha`, it computes

`D^alpha(w F) = sum_{beta<=alpha} binom(alpha,beta) D^beta w D^(alpha-beta) F`

with the exact Cartesian multi-index Leibniz coefficients. No finite
differences, fitting, sampled derivative reconstruction, or caller-supplied
final weighted derivative is used in this multiplication step.

The adapter fails closed unless stage, admitted scale `a_j`, exact evaluation
`q`, derivative order, `(A,B,p)` admission evidence, raw value provenance, and
cutoff-evaluator provenance all agree. It additionally requires the cutoff
jet's zero-order value to equal the previously evaluated cutoff weight and the
machine-derived zero-order weighted `(A,B,p)` values to equal the independent
finite-prefix contribution exactly.

## Truth boundary

The analytic product rule is machine-executed, but its two derivative tables
are still provider inputs. In particular this increment does **not** establish
that:

- the unweighted jet is generated from the actual Section 7/8 correction field;
- the cutoff jet is generated from the manuscript's implicit similarity
  coordinate and fixed cutoff;
- the infinite locally finite Eq. (9.21) field is constructed;
- the `t=1` endpoint is covered;
- a genuine Navier--Stokes residual artifact or Eq. (9.18) majorants exist; or
- the final smooth compact forcing, finite-energy certificate, or blow-up
  closure has been obtained.

Accordingly the downstream `Section9FinitePrefixJetCertificate` continues to
keep `analytic_derivatives_machine_derived_from_actual_corrections=false`,
`paper_fixed_cutoff_derivatives_machine_verified=false`,
`residual_artifact_ready=false`, and `paper_exact_velocity_available=false`.

## Next paper-exact dependency

The next substantive dependency is a real Section 7/8 correction-field exporter
that produces the unweighted `(A_j,B_j,p_j)` spacetime jets under the same
admission identity, together with a manuscript-derived analytic jet for
`chi(a_j q(t,z))`. Once those inputs are genuine, this adapter removes the need
to trust a provider-supplied final weighted derivative table; the resulting
prefix jet can then be used by a later local-velocity/Navier--Stokes residual
constructor without weakening the provenance boundary.
