# Section 10 paper time-switch endpoint-jet transparency

Status: **analytic adapter only**. This increment is bound to the existing pinned `openai/NavierStokesAndEuler` `SmoothCutoffs.lean` source and does not construct a new time window.

The pinned theorem-facing source already certifies that `timeSwitch(1)=1` and that every positive-order time derivative of `timeSwitch` vanishes at `t=1` because `1 > 3/4`. For a supplied finite spacetime jet of a scalar or three-vector field `F`, the adapter now evaluates the time-direction Leibniz rule

`D_t^n D_x^a D_y^b D_z^c (timeSwitch * F) = sum_{k=0}^n binom(n,k) timeSwitch^(k) D_t^(n-k) D_x^a D_y^b D_z^c F`.

At `t=1` only the `k=0` term remains, so the official Section 10 time switch preserves the supplied endpoint jet exactly. Regression coverage includes a test-side hand expansion of `D_t^2 D_x(timeSwitch*F)` and vector componentwise preservation; incomplete and non-finite supplied jets fail closed.

This narrows the endpoint blocker but does not solve it. In particular the repository still has no actual Section 9 all-order one-sided endpoint jet, no proof of convergence of those jets, no smooth extension through `t=1`, and no endpoint residual closure. Therefore all of the following remain false:

- `section9_endpoint_jet_supplied_by_actual_construction`
- `section9_field_smooth_extension_through_t1_constructed`
- `endpoint_residual_closure_verified`
- `paper_exact_velocity_available`

The next substantive step is to obtain the one-sided Section 9 field/residual jets from the actual correction construction and prove the all-order endpoint compatibility needed for a smooth extension. The Section 10 time switch can no longer be treated as a mechanism that supplies that missing regularity: at the endpoint it is analytically transparent.
