# Hierarchy-owned PositiveAxis base provenance

## Scope

This increment closes one ownership seam in the Section 5 Eq. (5.7) coefficient matrices. `A0` and `A1` consume an order-zero `PositiveAxisBaseJet`; previously that structured base jet could be maintained independently by a caller. The new bridge projects the required rows only from coefficient order zero of `Section5LowerHistorySixthMixedHierarchy`.

## Paper location and identity

The displayed PositiveAxis system in Section 5, Eq. (5.7), uses the leading profiles in the matrix coefficients. The landed `positive_axis_A0` and `positive_axis_A1` implementation consumes exactly these base quantities:

- `phi_0`, `partial_X phi_0`, `partial_X^2 phi_0`, `partial_eta phi_0`;
- `U_0`, `partial_X U_0`, `partial_X^2 U_0`, `partial_eta U_0`;
- `beta_0 = V_0 / X`.

The strong Stage-2 hierarchy already owns fifth-mixed `phi_0`, sixth-mixed `U_0`, and an analytic Eq. (5.2) fifth-mixed `beta_0`. `hierarchy_owned_positive_axis_base_jet(...)` therefore performs only exact downward projection of `phi_0/U_0` and takes `beta_0.value` from the hierarchy's Eq. (5.2) adapter. `hierarchy_owned_positive_axis_matrices(...)` then delegates the matrix values to the landed exact `positive_axis_A0/A1` formulas.

## Ownership chain

The matrix-value path is now

`order-zero strong hierarchy -> phi_0/U_0 high mixed jets -> analytic Eq. (5.2) beta_0 -> PositiveAxisBaseJet -> A0/A1`.

No independent `PositiveAxisBaseJet`, `beta` table, sampled `V/X` division, finite difference, fit, generic cutoff, or caller repair jet enters the public hierarchy bridge.

The order-zero strong data themselves remain Issue #1/upstream inputs. This increment does not claim that those inputs have already been materialized from the paper's final Theorem 4.6 profile.

## Regression contract

Focused tests verify that the bridge is exactly the order-zero hierarchy projection, that `beta_0` is the hierarchy-owned Eq. (5.2) value, that the produced matrices are byte-for-byte equal to the landed `positive_axis_A0/A1` formulas evaluated on that projection, and that unrelated hierarchy types fail closed. Test fixtures are analytic regression data only, not paper coefficient data.

## Truth boundary

Status remains `formal-structure`.

This increment owns only the **value rows** of `A0/A1`. Differentiating `K W_n^(0)` through `eta^2` still requires higher `eta`/mixed derivatives of the order-zero base quantities; those must come from the same genuine leading-profile hierarchy and may not be supplied as an unrelated derivative table. This increment also does not establish `W_n^(1)`, the Picard sum/convergence, a finalized positive-order `(phi_n,U_n,Pi_n,V_n)`, arbitrary-order coefficient assembly, uniform `C[j,m]` bounds, recursive cutoff convergence, Proposition 5.3, full reconstruction, or paper-exact velocity.