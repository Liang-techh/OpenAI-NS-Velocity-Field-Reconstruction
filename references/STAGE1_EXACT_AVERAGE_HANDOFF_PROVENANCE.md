# Stage-1 exact radial-average handoff provenance

## Scope

This is one downstream-only Issue #1 / Stage-1 profile-assembly increment built from exact repository `main` `3cc438f29524acb068e40bf3a99a8293c9222628`.

It does not implement `AxisCoefficientSpace`, coefficient operators, `naturalRemainder`, `x2`, a fixed point, or any missing scaled field. Agent 6 remains the owner of that upstream backend lane. No toy coefficient, sampled fit, default-zero state, or numerical near-zero identity is introduced here.

The repository already exposes `NaturalProfileAssembly` with six supplied fixed-point/profile objects: scaled `phi`, `u`, `du_deta`, the separate scaled radial-average field `average`, scaled pressure, and the axis pressure. Its Lean-faithful `radial_average` reconstruction is

`U_*(eta) + Lambda^-1 * B(Lambda X, eta)`.

Before this increment, `NaturalProfileAssembly.to_leading_profile()` discarded that already supplied `B` identity at the adapter boundary. `LeadingProfile.radial_average_U()` then fell back to a 32-point numerical Gauss-Legendre integral of `U`. That fallback is useful for arbitrary caller-supplied profiles, but it is not the correct handoff once the genuine natural fixed-point backend supplies the separate average field proved by the scaled system.

This increment therefore wires `NaturalProfileAssembly.radial_average` directly into `LeadingProfile.average_U`. The generic quadrature remains available for unrelated profiles, while the natural-profile path preserves the backend-owned average field exactly at the profile-interface level.

## Pinned formal source

Pinned source revision: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The relevant `NavierStokes/NaturalProfile.lean` structure is:

- `affineProfile`, which reconstructs the unscaled average by the same affine rescaling used for the axial field;
- `IsNaturalSolution.average_equation` and `IsNaturalSolution.average_integral`, where the true average is a separate solution field satisfying `X V(X,eta) = integral_0^X U(s,eta) ds`;
- `transportW_reconstruct`, which explicitly consumes the scaled average field `B` through `affineProfile (NaturalAxisData.U j) (1/Lambda) Lambda B`.

The adapter must therefore preserve `B`; silently replacing it by a fresh quadrature approximation is unnecessary identity loss at the downstream interface.

## Verification

`tests/test_natural_axis.py` now checks both the direct rescaling identity and a fail-fast regression in which the scaled axial `u` callable raises if sampled. `LeadingProfile.radial_average_U(...)` must still succeed from the supplied scaled `average` field, proving that the natural-profile handoff does not fall back to quadrature over `U`.

This verification is identity plumbing only. It does not certify that caller-supplied regression functions are the paper fixed point.

## Remaining boundary

The eta derivative of the average is intentionally not manufactured. `LeadingProfile.average_dU_deta` remains unresolved until a genuine backend/theorem handoff supplies the corresponding derivative field or exact certificate.

Stage 1 remains `formal-structure`. Complete `AxisCoefficientSpace` closure, genuine `naturalRemainder(x1)`, `x2`, the fixed point, actual `phi/u/average/pressure` materialization, `NaturalProfileAssembly` instantiation on that fixed point, regular-inner-core / heat-exterior matching, and support/moment/matching/cone closure remain open.

`full_reconstruction=false` and `paper_exact_velocity_available=false` remain mandatory.
