# Eqs. (5.7)-(5.8) inner-solver provenance

Status: **formal-structure**.

## Paper mapping

Section 5, Lemma 5.1 rewrites the positive-order inner coefficient equations in the radial variable `xi = sqrt(X)` as

`d_xi W_n + xi^-1 diag(0,0,2,0,3,1) W_n = A0 W_n + A1 d_eta W_n + f_n`, with `W_n(0,eta)=0`.

The manuscript then inverts the singular diagonal operator with

`(G g)_i(xi) = integral_0^xi (s/xi)^c_i g_i(s) ds`, `c=(0,0,2,0,3,1)`,

and writes the Picard operator as `K = G(A0 + A1 d_eta)`. The solution is represented by the convergent series `W_n = sum_{k>=0} K^k G f_n` on the fixed inner interval.

The sparse block structure of `A1` implies that a nonzero composition of `k` factors contains at most

`p_k = ceil(k/2)`

parameter derivatives. On a smaller complex strip `S_{rho'}` with `Delta=rho-rho'>0`, the paper's Eq. (5.8) gives

`||K^k G f_n||_{S_{rho'}} <= C_n^(k+1) a^(k+1)/(k+1)! * max(1,p_k/Delta)^p_k`.

The factorial simplex-volume factor dominates the Cauchy derivative loss, so the manuscript observes that the `k`-th root of the right-hand side tends to zero and hence the Picard series converges for every finite inner radius `a`.

## Executable artifacts

`src/openai_ns_reconstruction/background_inner_solver.py` implements:

- the exact six singular exponents from Eq. (5.7);
- the radial inverse `G`, including the zero-axis value without numerical division by `xi`;
- one Picard-map application `G[A0 W + A1 d_eta W + f_n]`; and
- the first positive-order Picard term `G f_n`.

`src/openai_ns_reconstruction/background_positive_axis.py` now removes the earlier arbitrary-matrix gap at the coefficient-algebra layer. It directly transcribes the pinned `NavierStokes/PositiveAxisSystem.lean` formulas for the real `A0`, sparse `A1`, and `f_n` forcing from structured leading/base jets and strictly lower-order source jets, fixes `lambda_n=2nh`, and exposes callbacks compatible with `background_inner_solver.py`. Its `first_positive_order_term_from_jets` therefore evaluates `G f_n` for the displayed theorem forcing rather than for an arbitrary six-vector callback.

`src/openai_ns_reconstruction/background_picard_bounds.py` implements the independent analytic-control side of Lemma 5.1:

- the exact derivative count `p_k=ceil(k/2)`;
- the displayed Eq. (5.8) term majorant, evaluated stably in logarithms;
- a fail-closed complete-tail enclosure. Once `p_k >= max(1,Delta)`, the code uses `(1+1/p)^p < 3` and `p+1 <= (k+3)/2` to obtain the conservative two-step bound
  `B_{k+2}/B_k <= 3(C_n a)^2/[2 Delta (k+2)]`;
- even/odd geometric summation then gives `sum_{j>=N} B_j <= (B_N+B_{N+1})/(1-r_N)` whenever the certified two-step ratio `r_N<1`; and
- a deterministic search for the first Picard truncation order whose analytic tail bound meets a requested tolerance. The order is chosen from Eq. (5.8), not fitted to sampled residuals.

`tests/test_background_inner_solver.py` cross-checks `G` against independently integrated polynomial primitives and the generic solver plumbing.

`tests/test_background_positive_axis.py` independently reconstructs the expanded profile right-hand sides corresponding to `PositiveAxisSystem.matrixRHS_jet`, checks the pinned `A1` sparsity, and verifies an actual `n=1` first term against closed-form integrals of the exact structured forcing.

`tests/test_background_picard_bounds.py` checks the displayed Eq. (5.8) arithmetic, compares the majorant with an exactly solvable scalar Picard specialization `K^kGf=lambda^k xi^(k+1)/(k+1)!`, independently sums a long finite segment of the paper majorants against the two-step tail enclosure, and verifies fail-closed strip/ratio/search-cap behavior.

## Boundary / blocker

The exact **formulas** for `A0`, `A1`, and `f_n` are now executable, but their structured `BaseJet`/`SourceJet` providers are not yet the paper-exact hierarchy. They must be derived from the materialized Section-4 leading profile and finalized lower orders, including the genuine lower convolutions and regular previous `Omega/X`. Issue #1 still does not provide that complete materialized leading profile. Caller-supplied jets therefore must not be labeled paper-exact.

The Eq. (5.8) module remains an analytic implication from independently certified inputs; it is not a proof that a supplied `C_n`, `rho`, or `rho'` is valid for the manuscript hierarchy. Its runtime evaluations use binary64 `log`/`lgamma`, and the singular inverse uses finite Gauss-Legendre quadrature, so these are reproducible numerical evaluations rather than interval-arithmetic or Lean proof objects.

Stage 2 therefore remains **formal-structure**, and `paper_exact_velocity_available` must remain `false`. The next profile-dependent step is to populate the landed positive-axis matrix/source adapter from the actual leading/lower-order recursion, certify its common analytic strip and `C_n`, and execute the converged Eq. (5.7) hierarchy before Lemma 5.2 compact repair.
