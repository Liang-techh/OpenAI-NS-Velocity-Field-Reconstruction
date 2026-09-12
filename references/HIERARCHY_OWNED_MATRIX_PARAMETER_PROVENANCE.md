# Hierarchy-owned Eq. (5.7) matrix parameter jets

## Scope

This increment analytically differentiates the displayed Section 5 PositiveAxis matrices `A0` and `A1` with respect to the similarity parameter `eta`, at fixed Picard radius `xi` (`X = xi^2`). It is a prerequisite for differentiating the first nontrivial hierarchy-owned Picard iterate landed in PR #184.

The production entry point is `hierarchy_owned_positive_axis_matrix_parameter_jets(...)` in `src/openai_ns_reconstruction/background_repaired_history_matrix_parameter.py`.

## Paper / formal source

The matrix values remain the already-landed transcription of Eq. (5.7) / the pinned PositiveAxis system in `background_positive_axis.py`, corresponding to the repository's official Lean modules:

- `NavierStokes/PositiveAxisSystem.lean`
- `NavierStokes/PositiveAxisExistence.lean`

The derivative uses only the order-zero hierarchy-owned second `(X, eta)` jets of `phi_0`, `U_0`, and the regular flux `beta_0 = V_0/X`. In particular `beta_0` and `partial_eta beta_0` are still generated through the analytic Eq. (5.2) regular-flux adapter; no independent beta derivative table is accepted.

## Analytic differentiation

At fixed `xi`, the implementation differentiates the exact scalar factors appearing in `A0`/`A1`, including

- `ell = 1 - 2 h eta^2`,
- `edge = 1 - eta^2`,
- `M = 1 - 2 eta U_0`,
- `R = M/ell + beta_0`,
- `Gphi = X partial_X phi_0 + phi_0`,
- `Gu = X partial_X U_0`, and
- the displayed `axialValue` quotient.

No production finite difference or interpolation is used.

## Regression

`tests/test_background_repaired_history_matrix_parameter.py` compares both analytic matrix derivatives against centered-`eta` finite differences of the pre-existing landed `hierarchy.positive_axis_fields(order).A0/A1` value path, at the axis and at positive `xi`. The finite difference exists only as an independent test oracle. The matrix values themselves are required to match the existing landed formulas exactly.

The test also verifies fail-closed behavior for an unowned recursive order and for a non-hierarchy input.

## Truth boundary

This is Stage-2 `formal-structure` infrastructure only. It does **not** establish `partial_eta W_n^(1)`, Picard convergence, a finalized positive-order coefficient, the recursive cutoff-scale schedule, Proposition 5.3 all-jets residual decay, or paper-exact velocity.

The reason `partial_eta W_n^(1)` is not claimed is structural: differentiating

`W_n^(1) = G(A0 W_n^(0) + A1 partial_eta W_n^(0) + f_n)`

also requires `partial_eta^2 W_n^(0)` through the term `A1 partial_eta^2 W_n^(0)`. The repository currently owns only the first `eta` jet of `W_n^(0)`. A hierarchy-owned second `eta` jet of the forcing / first Picard term is therefore the next derivative debt.

The genuine Issue #1 leading mixed profile remains upstream. Caller-supplied coefficients, generic cutoffs, sampled derivatives, and finite-order numerical fits are not promoted to paper-exact data. `paper_exact_velocity_available=false` remains mandatory.
