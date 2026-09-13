# Actual-schedule conditional profile-budget provenance

Status: **formal-structure only; conditional analytic budget**.

Artifacts:

- `src/openai_ns_reconstruction/axis_coefficient_profile_budget.py`
- `references/provenance_manifest_addendum_axis_coefficient_profile_budget.json`

The factory `actual_schedule_profile_budget(solver)` accepts only the genuine
`FormalAxisCoefficientSolverState`.  It reuses
`solver.x1.reference.reference.data` and `j` to rerun
`diagnose_actual_schedule_scale_chain_wide`.  It fails closed on an upstream
obstruction or on missing analytic norms, wide remainder, or wide scale.
Caller-fitted norm, remainder, Lambda, C, epsilon, and component bounds are not
accepted.

## Derived budget

The recomputed wide scale and symbolic C exponent must agree with the actual
amplitude scale carried by `solver.x1`; the recomputed Lambda must also agree
with `solver.Lambda`.  The epsilon taken from the actual analytic neighborhood
must agree with the source amplitude, the reference state, and the solver.

The wide remainder's `reference_norm_upper` is the existing actual-schedule
reference-pair max majorant, denoted `B0`.  The pinned contraction gate is
replayed with `NaturalPicardContractionCertificate.from_scale(scale)`.  If
`d` is its upward-rounded `one_step_radius_upper`, the stored profile budget is

```text
profile_norm_upper = ROUND_CEILING_96(B0 + d).
```

This is the common pair max majorant used for both angular and axial fields;
no separate fitted component norms are introduced.  The method
`conditional_tail_bound(...)` delegates this budget and the actual epsilon to
`axis_coefficient_radial_tail_bound`, which applies the pinned coefficient
weight and its derived negative-binomial radial tail estimate.

## Final ratio and pressure budgets

The same actual analytic chain exposes
`diagnostic.narrow_prefix.analytic_norms.amplitude_norm_upper`.  This is the
canonical normalized-amplitude bound from the componentwise Cauchy certificate,
not a caller-supplied constant.  The implementation imports that actual bound
into Decimal without binary64 under-rounding and forms every product in a
96-digit `ROUND_CEILING` context:

```text
M_a = actual analytic amplitude_norm_upper
M_phi = profile_norm_upper
M_F = ROUND_CEILING_96(64 * M_a * M_phi)
M_P = ROUND_CEILING_96(80 * 64 * M_F^2).
```

The factors are pinned by `AxisOperators.lean:326-331`, where
`norm_product_le <= 64`, and `AxisOperators.lean:467-483`, where
`norm_primitive_le <= 80`.  The assumptions are `epsilon > 0`, a common
`AxisSpace` coefficient norm for `phi`, and the actual analytic amplitude norm.
The pressure budget is for the scaled increment
`P = primitive(F^2)`; an independent axis pressure datum `P0` is outside this
coefficient budget.

`conditional_angular_ratio_tail_bound(...)` delegates the exact radial-tail
majorant with `M_F`, and `conditional_pressure_tail_bound(...)` delegates it
with `M_P`.  The physical factor `E = sqrt(2 * X) * F` is deliberately outside
this generic coefficient-weight budget and must be handled by profile assembly.

## Formal source and assumptions

The source pin is `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/AxisContraction.lean:347-349` defines the reference pair
  `x0 = (S one, -1/2 * J1(inverseL * zStar))`;
- `NavierStokes/AxisContraction.lean:502-511` supplies the fixed-point distance
  consequence used by the one-step radius;
- `NavierStokes/AxisWeightEstimates.lean:24-27` and
  `NavierStokes/AxisCoefficientSpace.lean:400-414` provide the coefficient
  weight and conditional norm-to-coefficient implication consumed by the
  delegated radial-tail module;
- `NavierStokes/AxisOperators.lean:326-331` gives the product norm bound 64;
- `NavierStokes/AxisOperators.lean:467-483` gives the primitive norm bound 80.

The existing `stage1_scale_chain_wide` and
`NaturalPicardContractionCertificate` are executable formal-structure numeric
majorants.  Their inputs still require the validity of the analytic/operator
norm inequalities, and the budget additionally requires identifying the
implemented compatible coefficient jets with the theorem's `AxisSpace`
fixed-point object.  Those assumptions are not proved by this factory.

## Truth and numerical boundary

The bounded implementation check was `python -m py_compile
src/openai_ns_reconstruction/axis_coefficient_profile_budget.py`.  A single
actual-schedule factory smoke constructed the formal solver and budget and
produced finite Decimal96 bounds for `M_a`, `M_F`, `M_P`, and both delegated
tails.  No focused pytest target was added or run for this extension.

The returned state reports `conditional_on_axis_space_identification=True`,
`global_axis_norm_certified=False`, `includes_coefficient_roundoff=False`, and
`paper_exact=False`.  The budget therefore does not certify global
`AxisSpace` membership, the formal solver's fixed point, convergence, a
paper-exact velocity, or a physical profile.  Decimal96 upward rounding of
the scalar budgets does not cover roundoff in evaluated coefficient jets, and
the delegated radial-tail bounds remain conditional on the supplied global
norm identification.

The repository base for this addendum is main commit
`67ab5ecd32f80225df8bce2760b943c20a7d3d7b`.
