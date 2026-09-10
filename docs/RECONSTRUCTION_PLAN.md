# Reconstruction plan and truth-status ledger

Updated 2026-09-10. The target remains an executable counterpart of the **entire published
velocity construction**, not merely the blow-up exponent. Every construction layer needs
source-mapped choices and independent verification. The goal has not been downgraded.

## Status convention

- **paper-exact:** direct published formula/algorithm; floating evaluation is still numerical.
- **formal-structure:** the structure is encoded but some profiles, parameters or recursive choices are missing.
- **diagnostic-only:** a numerical checker or experiment, not a proof construction.
- **toy:** illustrative data that must never be called the OpenAI counterexample.
- **pending:** no runtime implementation for the claimed layer.

`references/provenance_manifest.json` preserves the stable source/layer IDs introduced by
the parallel provenance work. The installed CLI additionally exposes stages 0–8 separately.
Tests check source-pin consistency. Neither metadata nor a caller's `paper_exact` flag proves
correctness. A green test run never automatically upgrades a stage.

## Stage 0 — similarity geometry: implemented formula

`coordinates.py` implements Eq. (4.1):

```text
A=1/2+h, D=1/2-h, tau=1-t
z=q^D eta, tau=q(1-eta^2), X=r^2/(2q)
q-z^2 q^(2h)=tau
```

The root solve is dimensionless and relative-scaled. `solve_q_from_tau` and
`similarity_coordinates_from_tau` avoid `1-tau` rounding; `d=tau/q` avoids subtracting
nearly equal numbers. Tests cover manufactured roots down to q=1e-200.
This is not a guarantee for every representable input or arbitrary precision.

## Stage 1 — leading velocity/profile: partial

`profiles.py` and `velocity.py` implement Eqs. (4.3)–(4.7):

```text
u_theta=q^(-A) E, u_z=q^(-A) U, r u_r=V0
AU=(1/X) integral_0^X U(x,eta) dx
V0=X/L [2 eta U-2D eta AU-(1-eta^2) d_eta AU]
L=1-2h eta^2, Pi_X=E^2/(2X)
```

The radial averages now use cached Gauss–Legendre nodes or user-supplied exact averages.
Regular-axis profile data can provide `F` with `E=sqrt(2X)F`; invalid nonzero axis swirl
is rejected instead of silently erased.

**New constructed component:** `heat_exterior.py` implements Appendix A.6, Eqs. (A.32)–(A.38):
H and its derivatives, exterior K(r,t), E_heat(X,eta), and centrifugal pressure normalized
at infinity. It includes derivative/ODE checks and independently tested NS balance for
this exterior-only flow. This component is not smooth at r=0 and cannot stand in for the core.

**Still required:** translate Theorem 4.6 and Appendices A/B/C into a regular inner profile,
matching moments, compact moment corrections, shear modification and admissible cone
conditions. Emit all constants/choices and generated coefficient hashes. Do not splice a
Gaussian into the missing core and mark the result complete.

## Stage 2 — all-order background: assembly only

`background.py` implements lambda_n=2nh, streamfunction/potential assembly and cutoff sums.
The generic cutoff in `cutoffs.py` is C-infinity mathematically, replacing the old C1
smoothstep, but **is not the paper-selected recursive cutoff schedule**.

Still required: solve Eqs. (5.2)–(5.6) for each coefficient, preserve compact moments,
and select cutoffs from the paper's estimates.
Acceptance: reproducible coefficients plus increasing derivative-order evidence for the
predicted super-algebraic residual decay. Finite-depth plots alone are insufficient.

## Stage 3 — dyadic charts and transported waves: geometry only

`charts.py` implements Eqs. (6.1)–(6.6): Q=2^-ell, epsilon=Q^h, S*=ell^2,
chart/physical conversion, the fixed integer covering matrix and its powers, and normalized
velocity/pressure/residual scaling. Q is fixed per chart; do not differentiate it as q(z,t).
Integer matrices use Python integers rather than overflowing int64. Floating fast-phase
accuracy is not certified by exact integer matrix arithmetic.

Still required: slow labels, separated auxiliary-torus supports, transported
wavevectors/polarizations, phase dynamics and rounded frequencies. Arbitrary frequency
lists are not acceptable for `paper-exact` status.

## Stage 4 — oscillatory stress realization: pending

Implement Section 7: admissible cone, positive covariance decomposition, primary waves,
amplitude solve, exact divergence-free realization by vector potentials, curl remainder,
physical evaluation and support separation.
Acceptance: independently computed averaged quadratic flux realizes the requested stress
to the stated order, with support and divergence identities checked.

## Stage 5 — compact mean corrections: pending

Implement Section 8 radial inverses, moment-preserving compact corrections and the
five-equation defect solve. Preserve all moments and supports needed by the iteration.
Acceptance: the zero-auxiliary-average defect is cancelled without breaking those constraints.

## Stage 6 — residual-improvement iteration: pending

Implement Section 9, including the recursive correction sequence and summation (9.21):

```text
Delta u_j=curl(A_j)+B_j
A=A0+sum_j chi(a_j q) A_j
B e_theta=B0 e_theta+sum_j chi(a_j q) B_j
u_loc=curl(A)+B e_theta
```

Choose a_j from estimates, not plot fitting. Acceptance: residuals and every tested derivative
vanish to increasing order with stable truncation studies. Full convergence still needs
analytic or formal bounds, not only a numerical experiment.

## Stage 7 — final compact field and force: composition only

`local_field.py` implements u=curl(cA)+cB e_theta, Eq. (10.4). Apply the cutoff before curl.
The direct swirl and its localized product must be axisymmetric to preserve divergence.
Factories encode that restriction; tests include a non-axisymmetric counterexample.

Still required: the completed local field, paper-specific support/cutoffs, force
f=u_t+(u.grad)u-Delta u+grad p at viscosity one, and the Section 10 smooth extension through t=1.
A force defined numerically as the residual has **not** thereby been shown smooth.

## Stage 8 — independent verification: partial diagnostics

Implemented: finite-input and stencil checks, bounded-time differentiation, manufactured
solutions with independently specified force, space/time refinement, exterior heat balance,
leading singular-path slope checks, deterministic CSV/JSON output, hashes and an audit CLI.

Still required: interval or symbolic/AD backends, kinetic energy/compact support checks,
full corrected-field convergence, trajectories/visualizations of the actual field, and
source-mapped Lean theorem cross-checks. Upstream is pinned, but no Lean build is claimed.

## Next independently reviewable work packages

1. Complete the regular inner/outer leading-profile constructor and matching/cone tests.
2. Implement one real Section 5 recursive coefficient and its residual/moment tests before generalizing.
3. Extend fixed charts to transported phases and verified support separation.
4. Implement the Section 7/8 stress and mean correction operators on independently generated defects.
5. Add analytic/interval error certificates, then integrate the Section 9/10 correction sequence.

These are uncompleted work packages, not claims that background jobs have been started.
Avoid parallel edits to the same path; re-read main and preserve concurrent work before committing.
