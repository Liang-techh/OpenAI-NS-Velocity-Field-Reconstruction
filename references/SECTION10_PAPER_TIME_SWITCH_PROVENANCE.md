# Section 10 paper time-switch provenance

Status: **theorem-source binding / support certificate**. This increment does not construct the missing one-sided Section 9 field extension through `t = 1`, does not evaluate Mathlib's noncomputable bump numerically, and does not promote paper-exact velocity or forcing.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- File: `NavierStokes/SmoothCutoffs.lean`
- Definition: `NavierStokes.SmoothCutoffs.timeSwitch`

The pinned source defines

`timeSwitch(t) = 1 - scaledCutoff(4/3)(t)`

from the fixed Mathlib `ContDiffBump` cutoff. The Python reconstruction does **not** replace that definition by the representative numerical bump in `cutoffs.py`.

## Exact theorem facts admitted

The new adapter binds only to the following pinned formal facts:

- `timeSwitch_contDiff`: the switch is `C^∞`;
- `timeSwitch_zero_of_abs_le`: `timeSwitch(t)=0` when `|t| <= 3/8`;
- `timeSwitch_one_of_three_quarters_le`: `timeSwitch(t)=1` when `t >= 3/4`;
- `timeSwitch_eventually_one`: for every `t > 3/4`, the switch is locally equal to one;
- `timeSwitch_iteratedDeriv_late`: every positive-order iterated derivative vanishes for `t > 3/4`;
- `timeSwitch_iteratedDeriv_support_nonneg`: on `t >= 0`, a nonzero positive derivative can occur only in the collar `[3/8, 3/4]`.

Consequently, at `t=1` the switch itself is theorem-certified to equal one and every positive-order switch derivative is theorem-certified to vanish. Thus the paper's switch is inert near the endpoint; an arbitrary user-selected time window is not an acceptable substitute.

## Landed gate

`src/openai_ns_reconstruction/section10_paper_time_switch.py` requires the exact repository, revision, file, theorem symbols, and rational geometry `3/8` and `3/4`. Boundary queries accept only exact `int`/`Fraction` inputs so binary floating-point rounding cannot change which formal theorem applies.

The adapter deliberately returns no transition value inside the unresolved collar and does not evaluate the Mathlib bump. Focused regression rejects source/revision/symbol drift, generic time-window geometry, and floating-point theorem-boundary queries.

## Explicit non-claims / next blocker

The facts above concern the multiplier `timeSwitch` only. Being identically one near `t=1` does **not** itself extend the Section 9 velocity/pressure field through `t=1`; multiplication by this switch leaves the endpoint behavior unchanged.

Therefore all of the following remain false:

- `actual_mathlib_bump_numerically_evaluated`
- `lean_theorems_machine_replayed_in_python`
- `section9_field_smooth_extension_through_t1_constructed`
- `endpoint_residual_closure_verified`
- `paper_exact_velocity_available`

The next truth-raising step on the time-localization path is the actual one-sided Section 9 field/residual extension through `t=1`, using the paper's all-order endpoint estimates. Final velocity and forcing must not be emitted before that construction and its independent closure checks exist.
