# Eq. (5.15) repaired-profile reconstruction provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: 50-52.

Implemented identities:

- Lemma 5.2, Eq. (5.14): after the compact `U_n`/`E_n` corrections, set `phi_n = C E_n/R` with the smooth axis extension inherited from the inner solution.
- Eq. (5.15):
  - `F_n(X,eta) = integral_0^X U_n(x,eta) dx`;
  - `V_n(X,eta) = - integral_0^X Z_{-A,n} U_n(x,eta) dx`;
  - `Pi_n(X,eta) = integral_0^X [C^-2 sum_{i+j=n} phi_i phi_j - Omega_{n-1}/(2x)] dx`.
- Eq. (5.2): the `V_n` integral is evaluated in the equivalent axis-regular form
  `V_n = [2 eta X U_n - 2 eta (D+lambda_n) F_n - d d_eta F_n]/L`, so production code never divides by `X` at the axis.
- Eq. (5.5): the pressure integral reuses the already-landed regular `Omega_{n-1}/X` pressure-row implementation.

## Executable mapping

`src/openai_ns_reconstruction/background_extension.py` provides `reconstruct_eq_5_15(...)` and the immutable `Eq515Reconstruction` result. The caller supplies the already-repaired analytic `U_n`, `d_eta U_n`, `phi_0,...,phi_n`, and the regular lower-order source `Omega_{n-1}/X`.

`tests/test_background_extension.py` checks the implementation in two independent ways:

- a polynomial positive-order profile is integrated directly with SciPy using the paper's unreduced `Z_b` operator from Eq. (4.2), rather than the production Eq. (5.2) rearrangement, and separately re-integrates the Eq. (5.15) pressure integrand;
- a second interoperability test feeds actual outputs of `Lemma52MomentRepair` into the reconstruction path and cross-checks `F_n` and `Pi_n` with adaptive quadrature.

The tests also verify the zero-axis forward-integration datum and fail closed on invalid orders, coordinates, normalization, missing `phi_n`, and nonfinite supplied profiles.

## Truth boundary

This remains **formal-structure / solver infrastructure**. The module does not construct the paper's order-`n` inner solution, does not derive `A0/A1/f_n`, does not generate the true eta-dependent Lemma-5.2 repair coefficients or their derivatives, and does not construct the upstream `Omega_{n-1}/X` hierarchy. Those data remain explicit dependencies.

Finite Gauss-Legendre quadrature is used to evaluate the forward integrals. Agreement with independent adaptive quadrature is numerical verification, not a uniform analytic support/stress certificate. In particular, this increment does not prove the Step-3 support closure in Lemma 5.2, the all-order coefficient bounds (5.17), the infinite cutoff/local-finiteness argument, or Proposition 5.3 residual decay.

No paper-exact gate is changed and `paper_exact_velocity_available` remains false.
