# Lemma 5.2 compact moment-repair provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Implemented formulas:

- Lemma 5.2, Step 1 and Eq. (5.14): add two compact `U_n` bumps and three compact `E_n` bumps in the reserved positive-order radial patch, fixed across orders.
- Eqs. (5.10)-(5.11): the five total radial moments that must vanish.
- Lemma 5.2, Step 2 and Eq. (5.16): with `E_0=e_* f(eta) R^(-1-2 lambda)` on the reserved patch, the five affine moment conditions split into `B_U alpha_n = -d_{U,n}` and `B_E beta_n = -d_{E,n}`.
- Lemma A.1: distinct exponents and ordered disjoint bump supports give invertibility of the two moment matrices.

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: 49-52.

## Executable mapping

`src/openai_ns_reconstruction/background_moment_repair.py` provides:

- an explicit nonnegative unit-mass `C^infinity` bump on each positive `R` interval;
- a fixed five-bump geometry with pairwise disjoint supports;
- the manuscript moment matrices with exponents `p=(1,1-2 lambda)` and `s=(2,-2-2 lambda,-2 lambda)`;
- the exact Eq. (5.16) block solve for `alpha_n(eta)` and `beta_n(eta)`;
- Eq. (5.14) value adapters for the compact corrections to `U_n` and `E_n`;
- fail-closed checks for malformed geometry, zero `e_* f(eta)`, nonfinite data, or a numerically singular matrix.

`tests/test_background_moment_repair.py` independently re-integrates the five correction moments with SciPy adaptive quadrature rather than the production Gauss-Legendre rule and checks that all five repaired moments vanish to floating-point tolerance.

## Truth boundary

This is a **formal-structure / solver-infrastructure** increment, not a paper-exact completed positive-order coefficient.

The unrepaired moment vector `m_n^0(eta)` and the patch factor `e_* f(eta)` are explicit inputs. They are not yet generated from the fully materialized paper profiles. Numerical matrix inversion is not a formal proof of Lemma A.1, and sampled calls do not establish smoothness or a uniform lower bound for `f` over `|eta|<=1`.

The module also does not yet reconstruct the full repaired `F_n,V_n,Pi_n` family from Eq. (5.15), prove the support/stress conclusions of Lemma 5.2, choose the recursive `q` cutoff scales, or prove Proposition 5.3 residual decay. `paper_exact_velocity_available` therefore remains false.
