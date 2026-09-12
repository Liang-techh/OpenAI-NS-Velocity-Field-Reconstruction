# Section 9 cutoff-composition jet provenance

## Scope

This increment narrows one Issue #4 derivative trust boundary for the positive-stage Eq. (9.21) weight `chi(a_j q)`.  The previously landed `section9_weighted_correction_jet.py` accepts a full spacetime derivative table for that scalar weight.  `section9_cutoff_composition_jet.py` instead accepts a finite spacetime jet of `q` and the one-variable values `chi^(k)(a_j q)`, then derives the full weight jet by normalized multivariate Taylor composition.

For normalized coefficients `f_alpha = D^alpha f / alpha!`, put `delta = a_j (q-q0)`.  The implementation truncates only terms above the requested total derivative order and computes

`chi(a_j q0 + delta) = sum_k chi^(k)(a_j q0) delta^k / k!`.

This is the finite multivariate chain/Faa-di-Bruno algebra.  No finite differences, fitted coefficients, derivative sampling, or reuse of the final weighted derivative table is involved.

## Source identities and cross-checks

The adapter is tied to an already-admitted Section 9 stage and an already-landed finite-prefix evaluation.  It requires:

- the q-jet zero-order value and exact rational identity to match the finite-prefix evaluation q;
- the scalar cutoff derivative sample to be evaluated at exactly `a_j q`;
- the scalar cutoff zero-order value to match the independently evaluated finite-prefix cutoff weight;
- the scalar cutoff provenance to match the finite-prefix cutoff evaluator provenance; and
- complete derivative coverage through the requested finite total order.

The regression test independently checks the literal mixed derivative

`D_tx chi(aq) = chi''(aq) (a q_t)(a q_x) + chi'(aq) a q_tx`

rather than recomputing it through the Taylor-composition helper.

## Truth boundary

This file does **not** establish the manuscript's actual cutoff-weight derivatives.  Both the spacetime q-jet and the one-variable `chi^(k)` data remain provider inputs.  In particular:

- `similarity_coordinate_jet_machine_derived_from_eq_4_1 = false`;
- `paper_fixed_cutoff_derivatives_machine_verified = false`; and
- `paper_exact_velocity_available = false`.

The adapter does not construct the infinite Eq. (9.21) field, the `t=1` endpoint extension, a genuine Navier-Stokes residual artifact, smooth compact forcing, finite-energy control, or blow-up closure.

The next substantive step is to replace the q-jet input by an all-order jet generated from the actual Eq. (4.1) similarity coordinate and to bind the one-variable cutoff derivatives to the manuscript's fixed cutoff source.  Until both happen, this module is analytic plumbing only.
