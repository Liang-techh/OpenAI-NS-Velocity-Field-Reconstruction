# Hierarchy-owned PositiveAxis forcing through third eta derivative

## Scope

This increment closes one finite derivative seam in Issue #2 / Section 5.  It composes the hierarchy-owned strict-lower source jet with the already-landed displayed PositiveAxisSystem forcing formula from Eq. (5.7).  It does **not** solve the positive-order fixed point, construct a repaired coefficient, or prove any all-order convergence statement.

Current integration baseline checked before this increment: `main` at `13f9c53bbe1a68ddd3fdaf1805e7f8b6bbd38b5d`.  The implementation is intentionally stacked on PR #319 head `d60b7b9e41ccdfe88880a287bee0702355a3fbd2`, because the new constructor consumes #319's hierarchy-owned `partial_eta^3 actualLowerSource` and #319 is still awaiting its authoritative Actions run.  No claim is made that either stacked increment is already integrated into `main`.

## Pinned formula

For positive order `n`, squared radius `X=xi^2`, and

`p_n = C^-2 source.pressure_product - source.omega_quotient/2`,

we use the displayed PositiveAxisSystem forcing

`f_n = (0, 0, 0, 2 xi p_n, 2 source.angular, 2 source.axial - 4 X (eta/ell) p_n)`

with `ell = 1 - 2 h eta^2`.

The value row is delegated exactly to the existing `positive_axis_forcing`.  For eta derivatives one through three, the only new scalar geometry is the analytic jet of

`q(eta)=eta/ell`:

- `q' = (1 + 2 h eta^2)/ell^2`,
- `q'' = 4 h eta (3 + 2 h eta^2)/ell^3`,
- `q''' = 12 h/ell^2 + 192 h^2 eta^2/ell^4`.

The product `q p_n` is differentiated with exact binomial Leibniz rows.  Because `p_n` is linear in the source pressure-product and Omega/X entries, its eta derivatives are obtained by applying the same landed `pressure_source` linear map to each hierarchy-owned source derivative row.

## Ownership and fail-closed boundary

Production accepts only `Section5LowerHistorySixthMixedHierarchy`.  It calls `hierarchy_owned_lower_source_third_parameter_jet(...)`; callers cannot supply an independent `SourceJet` derivative table or an unrelated normalization `C`.  `h` and `C` are taken from the same hierarchy that owns repaired fifth-mixed phi, repaired sixth-mixed U, and the Eq. (5.2)-derived beta chain.

No finite difference, fit, generic cutoff, sampled `V/X` division, or hand-filled forcing derivative appears in production.  Binary64 overflow/nonfinite output fails closed.

## Regression

The focused regression reuses the coherent analytic strong-jet fixture already used for #319.  It checks:

1. the forcing value row is exactly the landed `positive_axis_forcing` value;
2. the analytic third eta row agrees with a centered derivative of the independently evaluated analytic second eta forcing row;
3. the axis `xi=0` remains finite and the pressure radial component keeps its exact zero factor there; and
4. a non-strong hierarchy type is rejected.

The centered difference is a **test-only oracle**.  The analytic fixture is not paper coefficient data.

## Truth boundary / Agent 7 split

This remains Stage 2 `formal-structure` with `full_reconstruction=false` and `paper_exact_velocity_available=false`.  In particular this increment does not provide `W_n^(0)'''`, `W_n^(1)''`, a converged PositiveAxisSystem/Picard coefficient, a materialized repaired `(phi_n,U_n,Pi_n,V_n)`, an all-order hierarchy, Proposition 5.3, recursive-cutoff convergence, or residual all-jets-flatness.

Agent 7 retains ownership of finite/stagewise `C[j,m]`, recursive SlowBorel/DiagonalScale scheduling, common-support binding, target-order arithmetic, and residual-convergence certificates.  This increment touches none of those paths.
