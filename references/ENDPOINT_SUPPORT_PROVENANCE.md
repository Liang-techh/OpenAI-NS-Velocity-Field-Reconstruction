# Section 10 endpoint-support certificate provenance

Status: **formal-structure**. This increment does not certify the final paper velocity or force.

## Pinned source

Official formalization: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The relevant algebra is in `NavierStokes/SpatialBorelExtension.lean`:

- `extension_zero_of_coefficients_zero` proves that if every spatial Taylor--Borel coefficient vanishes at a fixed spatial point, then the extension vanishes there for every parameter value;
- `rightExtension` supplies the endpoint-shifted version used by the spacetime glue;
- the existing repository `DoublingEnvelope` is the executable local-finiteness mechanism for the future series.

`NavierStokes/SupportedParameterExtension.lean` uses the same theorem to show that zero spatial fibers persist through a Taylor--Borel gluing, making explicit that support preservation comes from zero endpoint coefficients rather than from enlarging the support set.

## Executable mapping

`src/openai_ns_reconstruction/endpoint_support.py` adds
`formal_force_support_point_certificate(...)`.  At a point outside the official Section 10 closed cylinder

`r^2 <= 1/16, |z| <= 1/4`,

it follows the same branch convention as
`traced_residual.formal_force_from_full_spacetime_jets`:

- for `t<T`, it checks the supplied past residual value directly;
- for `t=T`, only the degree-zero endpoint coefficient can contribute;
- for `t>T`, it uses the landed doubling scales and cutoff support to enumerate exactly the finite set of normal Taylor coefficients that can contribute at that time, requires exact zero for every such coefficient, and cross-checks that the executable Borel value is exactly zero.

No tolerance is used to decide support.  A malformed full spacetime tensor, invalid scale, or nonzero active coefficient fails closed.

## Regression boundary

`tests/test_endpoint_support.py` checks all three branches.  The future regression uses local scale `1`, hence the doubling envelope `1,2,4,8,...`; at `t-T=0.2`, only degrees `0,1,2` can contribute and the test deliberately raises if degree `3` or later is sampled.  A separate regression injects a nonzero active degree-one normal coefficient and confirms that the certificate does not pass.  Another test confirms that a nonzero past residual also fails closed.

These manufactured tensors are regression fixtures only.  They are not endpoint tensors from the paper construction.

## What remains before paper-exact support

This increment does **not** prove that the actual Section 9/10 localized Navier--Stokes residual or all of its endpoint normal jets vanish outside the support cylinder.  That requires the unresolved actual local field, its smooth localized residual, and the genuine locally-uniform endpoint derivative limits.  It also does not turn finitely many point queries into a global compact-support theorem.

Accordingly Stage 7 remains **formal-structure**, `full_reconstruction` remains `false`, and `paper_exact_velocity_available` must remain `false` until the actual upstream field, residual limits, smooth endpoint force, blow-up path, and energy obligations are closed.
