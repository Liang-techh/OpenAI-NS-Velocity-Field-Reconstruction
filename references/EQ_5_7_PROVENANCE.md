# Eq. (5.7) inner-solver provenance

Status: **formal-structure**.

## Paper mapping

Section 5, Lemma 5.1 rewrites the positive-order inner coefficient equations in the radial variable `xi = sqrt(X)` as

`d_xi W_n + xi^-1 diag(0,0,2,0,3,1) W_n = A0 W_n + A1 d_eta W_n + f_n`, with `W_n(0,eta)=0`.

The manuscript then inverts the singular diagonal operator with

`(G g)_i(xi) = integral_0^xi (s/xi)^c_i g_i(s) ds`, `c=(0,0,2,0,3,1)`,

and writes the Picard operator as `K = G(A0 + A1 d_eta)`.  The solution is represented by the convergent series `W_n = sum_{k>=0} K^k G f_n` on the fixed inner interval.

## Executable artifact

`src/openai_ns_reconstruction/background_inner_solver.py` implements:

- the exact six singular exponents from Eq. (5.7);
- the radial inverse `G`, including the zero-axis value without numerical division by `xi`;
- one Picard-map application `G[A0 W + A1 d_eta W + f_n]`;
- the genuine first positive-order Picard term `G f_n`.

`tests/test_background_inner_solver.py` cross-checks `G` against independently integrated polynomial primitives, verifies an actual `n=1` first Picard term, exercises a nonzero `K` contribution with an independent closed-form oracle, and checks fail-closed input validation.

## Boundary / blocker

This does **not** materialize the paper's profile-dependent matrices `A0`, `A1`, or lower-order source `f_n`.  Those depend on the leading profile and finalized lower-order coefficients.  Issue #1 still does not provide the paper-exact materialized leading profile, so caller-supplied matrices or sources must not be labeled paper-exact.  The new code is solver infrastructure for the Lemma-5.1 construction, not a completed positive-order coefficient and not evidence that Stage 2 is paper-exact.

The next Section-5 step is to derive the actual Eq. (5.3)-(5.6) coefficient matrices/source from materialized profile data, then iterate Eq. (5.7) to a certified solution on the common inner interval before applying Lemma 5.2 compact moment repair.
