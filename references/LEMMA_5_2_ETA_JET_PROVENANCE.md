# Lemma 5.2 eta-jet repair provenance

Status: **formal-structure**.

## Paper mapping

Section 5, Lemma 5.2 and Eqs. (5.14)-(5.16) repair each positive-order coefficient with two fixed compact `U_n` bumps and three fixed compact `E_n` bumps.  The already-landed `background_moment_repair.py` realizes the five constant moment matrices pointwise in `eta`.

For recursive use at the next order, the repair coefficients must vary smoothly with `eta`.  Since the bump geometry and the matrices `B_U`, `B_E` are independent of `eta`, all parameter dependence in Eq. (5.16) enters through the five unrepaired moments `m_i(eta)` and the patch factor

`p(eta) = e_* f(eta)`.

Writing ordinary derivatives rather than factorial-normalized Taylor coefficients, derivatives of the three ratios needed by Eq. (5.16) are propagated by the identity

`(m/p)^(n) = [m^(n) - sum_{j=1}^n binom(n,j) p^(j) (m/p)^(n-j)] / p`.

At each derivative order the same constant matrices are then solved:

- `B_U alpha^(n) = -[(m_0)^(n), (m_3/p)^(n)]`;
- `B_E beta^(n) = -[(m_1)^(n), (m_2/p)^(n)/2, -(m_4/p)^(n)]`.

The differentiated five corrected moments are recomputed with the Leibniz rule for the products `p B_U[1] alpha`, `p B_E[1] beta`, and `p B_E[2] beta`.

## Executable mapping

`src/openai_ns_reconstruction/background_moment_repair_jets.py` adds:

- a finite eta-jet adapter for supplied ordinary derivatives of the five unrepaired moments and `p=e_*f`;
- exact binomial-coefficient quotient and product derivative recurrences;
- derivative-order solves against the already constructed fixed `B_U` and `B_E` matrices;
- recomputation of all five differentiated corrected-moment identities;
- fail-closed checks for malformed/nonfinite jets, zero center patch factor, overflow, and a floating-point algebraic cancellation failure.

`tests/test_background_moment_repair_jets.py` does not manufacture the same quotient recurrence as its oracle.  It chooses independent polynomial `p(eta)`, `alpha(eta)`, and `beta(eta)`, generates the five unrepaired moment polynomials in the forward Eq. (5.16) direction using `numpy.polynomial.polynomial.polymul`, converts their Taylor coefficients to ordinary derivatives, and verifies that the production jet solver recovers the prescribed coefficient jets.  It also checks that order zero exactly bridges to the existing pointwise `Lemma52MomentRepair.solve` result and exercises fail-closed inputs.

## Truth boundary

This increment does **not** derive the supplied moment jets from the actual Eqs. (5.3)-(5.7) hierarchy.  It also does not prove that `e_* f(eta)` is uniformly separated from zero on an eta interval, does not infer analytic strip bounds from samples, and does not prove the full Lemma 5.2 support/stress conclusions required for the next recursive order.

The runtime matrix solves and cancellation checks are binary64 verification of finite algebraic jets, not Lean proof objects or interval certificates.  A green test therefore does not promote the result to paper-exact.

Stage 2 remains **formal-structure** and `paper_exact_velocity_available` remains `false`.  Once Issue #1 and the lower-order recursion provide genuine analytic moment/patch-factor jets, this adapter can carry the compact repair through eta derivatives without replacing the manuscript construction by pointwise fitting.
