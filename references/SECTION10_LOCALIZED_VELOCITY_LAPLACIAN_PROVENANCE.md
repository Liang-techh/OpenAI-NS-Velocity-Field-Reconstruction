# Section 9 -> Section 10 localized velocity Laplacian provenance

Status: **formal-structure analytic derivative adapter only**.

Base integration commit: `a5431dd5aac603fb0af6aa91f601d8ab40c8d66e`.

## Inputs already admitted on `main`

- `Section9FinitePrefixJetCertificate`: a complete provider-supplied finite spacetime jet for one admitted Eq. (9.21) prefix.  It is not the manuscript's constructed infinite correction sequence.
- `section9_section10_localized_velocity_spatial_jacobian.py`: analytic off-axis `grad u_cut` for
  `u_cut = c curl(A) + grad(c) x A + c B e_theta`.
- `spatial_localization.py`: the fixed Section 10 support/plateau geometry and analytic gradient/Hessian of the repository's explicit C-infinity representative.
- `section10_spatial_cutoff_viscous_adapter.py`: analytic `Delta c` and `grad(Delta c)` for that same executable representative.

The official geometry remains `r^2 < 1/32, |z| < 1/8` on the plateau and support inside `r^2 <= 1/16, |z| <= 1/4`.  Transition-collar values are **not** claimed pointwise equal to Mathlib's noncomputable `ContDiffBump`.

## New bounded identity

For `V = curl(A) + B e_theta`, `g = grad(c)`, `H = Hess(c)`, and `L = Delta c`, the adapter evaluates off-axis

`Delta u_cut = L V + 2 grad(V) g + c Delta V + grad(L) x A + 2 sum_j H[:,j] x partial_j A + g x Delta A`.

It requires total derivative order at least three because `Delta curl(A) = curl(Delta A)` uses third spatial derivatives of `A`.  `Delta(B e_theta)` is differentiated analytically in Cartesian coordinates for `r>0`; the symmetry axis fails closed pending the actual paper field's axis-regularity data.

Production code performs no finite differencing.  Regression validation finite-differences the complete existing `LocalizedField.velocity` as an independent oracle and separately hand-checks nonzero `Delta curl(A)` and `Delta(B e_theta)` manufactured values.

## Explicit non-claims

This adapter does **not** establish any of the following:

- actual Section 7/8 correction-field export;
- the infinite/locally-finite Eq. (9.21) construction;
- one-sided all-order convergence or smooth extension through `t=1`;
- a genuine Navier--Stokes residual artifact or smooth forcing;
- endpoint residual closure;
- the global symmetry-axis regularity needed by the final field;
- bounded-energy or blow-up closure;
- paper-exact velocity availability.

Accordingly `paper_exact_velocity_available=false`, `residual_artifact_ready=false`, and all actual/infinite/endpoint completion flags remain mandatory.  A future genuine residual constructor may consume this Laplacian only after the upstream correction field and endpoint construction are independently materialized.
