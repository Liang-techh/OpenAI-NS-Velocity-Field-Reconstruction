# Section 9 -> Section 10 convective acceleration provenance

Status: **formal-structure analytic operator adapter only**.

Base integration commit: `37a1fda0afe961014a55e38c6f2887846bf1b6f4`.

## Scope

This increment adds only the convective term required by a future genuine
localized Navier--Stokes residual.  It consumes the already-landed finite-prefix
Section 9 -> fixed Section 10 localized velocity and analytic off-axis spatial
Jacobian and evaluates

`(u_cut . grad) u_cut = (grad u_cut) u_cut`.

The spatial-Jacobian convention is fixed by the landed adapter: column `j` is
`partial_j u_cut`, so matrix-vector multiplication by `u_cut` implements
`sum_j u_j partial_j u_i` componentwise.  No caller-supplied cutoff, velocity,
or Jacobian is accepted and production code performs no finite differencing.

## Independent regression

The regression uses an analytic polynomial vector potential with zero swirl at
a point where both factors of the executable Section 10 cutoff are in their
transition collars.  As an independent oracle it constructs the complete
`LocalizedField` **without** an analytic cutoff gradient, so that the poloidal
piece is obtained by numerically curling `cA` as a whole.  It then freezes the
independently evaluated velocity `u(x0)` and centered-finite-differences the
complete localized velocity along the line `x0 + s u(x0)`.  This checks the
new convective acceleration without reusing the production spatial Jacobian.
The regression separately requires exact zero outside the official fixed
support and verifies the inherited second-order/truth gates.

## Truth boundary

The underlying `Section9FinitePrefixJetCertificate` is still a provider-supplied
finite prefix.  This increment does **not** establish:

- the actual Section 7/8 correction-field exporter;
- the infinite/locally-finite Eq. (9.21) field;
- one-sided all-order convergence or smooth extension through `t=1`;
- symmetry-axis regularity of the final paper field;
- a genuine Navier--Stokes residual artifact or smooth forcing;
- endpoint residual closure;
- bounded kinetic energy or blow-up closure;
- paper-exact velocity availability.

Accordingly `actual_section9_sequence_verified=false`,
`eq_9_21_infinite_sum_constructed=false`,
`section9_field_smooth_extension_through_t1_constructed=false`,
`residual_artifact_ready=false`, `endpoint_residual_closure_verified=false`, and
`paper_exact_velocity_available=false` remain mandatory.
