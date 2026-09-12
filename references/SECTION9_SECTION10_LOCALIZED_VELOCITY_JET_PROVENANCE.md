# Section 9 prefix jet -> Section 10 localized velocity provenance

Status: **finite-prefix analytic adapter only**. This increment connects two already-landed formal-structure paths and deliberately does not construct a Navier--Stokes residual or forcing.

## Inputs and source mapping

- `src/openai_ns_reconstruction/section9_eq921_prefix_jet.py` supplies an admitted finite Eq. (9.21) prefix jet `(A,B,p)` with complete physical spacetime multi-index coverage through a caller-selected finite order. Its own truth gate still records that the derivatives are not the actual infinite manuscript correction sequence.
- `src/openai_ns_reconstruction/local_field.py` fixes the local representation `u_loc = curl(A) + B e_theta` and the localized product rule `curl(cA)=c curl(A)+grad(c) x A`.
- `src/openai_ns_reconstruction/spatial_localization.py` fixes the executable Section 10 support/plateau representative and its matching analytic Cartesian gradient. The geometry is bound to the pinned `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd` `NavierStokes/SpatialLocalization.lean` support/plateau statements; transition-collar values remain this repository's explicit C-infinity representative rather than a claim of pointwise equality with Mathlib's noncomputable `ContDiffBump`.

`section9_section10_localized_velocity_jet.py` reads the first Cartesian derivatives of `A` directly from the prefix jet,

`curl(A) = (A_z,y - A_y,z, A_x,z - A_z,x, A_y,x - A_x,y)`,

forms `B e_theta`, and then evaluates

`u_cut = c curl(A) + grad(c) x A + c B e_theta`

with the fixed Section 10 `(c, grad(c))` pair. The production path does not finite-difference either `A` or the cutoff.

## Independent regression

The regression uses the manufactured polynomial `A=(0,0,xy+0.3x+0.2z)` and constant `B=0.3`. It checks the analytic curl by hand, then compares the adapter's localized velocity in the Section 10 transition collar against the pre-existing `LocalizedField` path with `analytic_curl=None`; that path finite-differences `curl(cA)` as a whole, so it does not reuse the analytic product-rule evaluation under test. A separate outside-support check confirms exact zero from the fixed cutoff/gradient pair.

## Truth boundary

The finite prefix certificate does not carry an independently verified spacetime base point, so the `(x,y,z,t)` binding in this adapter is caller-supplied metadata. More importantly, the actual Section 7/8 correction-field exporter, the infinite/locally-finite Eq. (9.21) sequence, the one-sided all-order `t->1-` field/residual limits, smooth extension through `t=1`, and genuine residual/forcing closure are still missing.

Accordingly the adapter keeps all of these false:

- `point_binding_to_actual_correction_field_verified`
- `actual_section9_sequence_verified`
- `eq_9_21_infinite_sum_constructed`
- `section9_field_smooth_extension_through_t1_constructed`
- `residual_artifact_ready`
- `endpoint_residual_closure_verified`
- `paper_exact_velocity_available`

This file is therefore a reusable analytic bridge toward the later genuine residual constructor, not evidence that the paper-exact velocity or forcing has been reconstructed.
