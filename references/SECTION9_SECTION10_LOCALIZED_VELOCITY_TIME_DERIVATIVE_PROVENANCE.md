# Section 9 -> Section 10 localized velocity time-derivative provenance

Status: **formal-structure**. This increment does not make the reconstructed velocity paper-exact.

## Mathematical source boundary

The adapter uses the already-landed finite-prefix Eq. (9.21) analytic jet together with the fixed Section 10 spatial localization algebra

`u_cut = curl(c A) + c B e_theta = c curl(A) + grad(c) x A + c B e_theta`.

Because the Section 10 spatial cutoff is time independent, differentiating this identity gives

`d_t u_cut = c curl(d_t A) + grad(c) x d_t A + c (d_t B) e_theta`.

The mixed derivatives needed by `curl(d_t A)` are read analytically from the prefix jet. No production finite difference is used.

## Repository increment

`src/openai_ns_reconstruction/section9_section10_localized_velocity_time_derivative.py` requires a complete finite-prefix spacetime jet through total order at least two, reuses the fixed Section 10 cutoff/gradient bridge, and returns the one-point analytic time derivative of the localized velocity.

`tests/test_section9_section10_localized_velocity_time_derivative.py` uses a manufactured time-dependent vector potential and an axisymmetric radial swirl amplitude. The regression independently finite-differences the *complete* localized velocity in time; the nested velocity path numerically curls `cA` as a whole and therefore does not reuse the production time-derivative product rule. A separate regression checks exact zero outside the official fixed spatial support.

## Explicit non-claims

The input jet is still provider supplied. This increment does not verify the actual Section 7/8 correction fields, construct the infinite/locally-finite Eq. (9.21) field, establish all-order one-sided convergence as `t -> 1-`, construct the smooth extension through `t=1`, build the genuine Navier--Stokes residual artifact, prove forcing smoothness/closure, or verify finite energy/blow-up.

The following therefore remain false:

- `actual_section9_sequence_verified`
- `eq_9_21_infinite_sum_constructed`
- `section9_field_smooth_extension_through_t1_constructed`
- `residual_artifact_ready`
- `endpoint_residual_closure_verified`
- `paper_exact_velocity_available`

This is one reusable analytic derivative adapter toward the future forcing calculation, not a forcing constructor.
