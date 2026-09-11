# Section 10 closed-past extension provenance

Status: **formal-structure**. This increment materializes the algebraic closed-past branch used immediately before the endpoint residual-limit theorem. It does not construct the missing paper-exact velocity/pressure, prove smoothness, or supply the all-order endpoint derivative limits required by `CandidateFromLimits`.

## Official source mapping

Pinned official Lean repository:

`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant modules and statements:

- `NavierStokes/SmoothCutoffs.lean`: constructs `timeSwitch`, with a zero germ around `t=0` and value one from `t>=3/4`.
- `NavierStokes/TimeLocalization.lean`: defines `activatedVelocity u = timeSwitch * u` and `activatedPressure p = timeSwitch * p`.
- `NavierStokes/PastExtension.lean`: defines `zeroBefore`, `pastVelocity`, `pastPressure`, and `pastResidual`; negative times are replaced by zero, while nonnegative times keep the activated fields. The zero germ at `t=0` is the theorem-side reason the branch is jointly smooth on the open past.
- `NavierStokes/CandidateFromLimits.lean`: consumes `PastExtension.pastResidual` when forming the traced residual whose actual full spacetime derivative limits are glued through `t=1`.

## Executable mapping

`src/openai_ns_reconstruction/past_extension.py` implements:

- `zero_before_vector` / `zero_before_scalar`: executable counterparts of `PastExtension.zeroBefore` in the repository's `(x,y,z,t)` calling convention;
- `past_velocity(u) = zero_before_vector(activated_velocity(u))`;
- `past_pressure(p) = zero_before_scalar(activated_pressure(p))`;
- `past_residual_numeric`: an independent finite-difference evaluation of the Navier--Stokes residual of those closed-past fields for `0<t<1`, with exact structural zero returned on `t<=0` and an open-endpoint temporal stencil near `t=1`.

The negative branch short-circuits the unresolved input functions. This matters because the official construction deliberately assumes no negative-time regularity of the original physical fields.

## Independent checks

`tests/test_past_extension.py` verifies:

- negative-time vector/scalar branches return zero without evaluating unavailable inputs;
- the combined past fields have the pinned zero germ through the executable `timeSwitch` plateau and recover the original fields after `t>=3/4`;
- the residual is exactly zero on `t<=0` without sampling unresolved inputs;
- on the positive transition collar, the residual of the closed-past fields agrees with the separately evaluated `TimeLocalization.activated_residual_formula_numeric` identity rather than with a copied implementation of `past_residual_numeric`;
- the diagnostic refuses evaluation at or beyond the singular endpoint.

## Truth boundary

This increment does **not** establish any of the following:

- equality of the executable time-switch transition values with Mathlib's noncomputable `ContDiffBump` pointwise in the collar;
- theorem-level `ContDiffOn` smoothness of `zeroBefore` or `pastResidual`;
- the full spacetime derivative recurrence of the actual residual;
- `TendstoLocallyUniformly` of every actual residual derivative as `t -> 1-`;
- the endpoint derivative majorants used by `SpatialBorelExtension`;
- a smooth global force through `t=1`;
- the unresolved paper-exact local velocity/pressure from earlier construction stages.

`past_residual_numeric` is diagnostic. In particular, defining a later force by this numerical residual and checking `R-f=0` would be tautological and is not treated as independent verification.

Therefore Stage 7 remains **formal-structure**, and `paper_exact_velocity_available` must remain `false`.
