# Axis coefficient parameter-primitive provenance

## Scope

This increment materializes one exact coefficient-space operator required by the
pinned natural-axis contraction: `AxisOperators.parameterPrimitive` on the
already landed actual `SchedulePressure`-derived coefficient eta-jets.

It is **not** a complete `AxisSpace` implementation and does not make the
leading profile paper-exact.

## Pinned upstream source

Repository: `openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant definitions:

- `NavierStokes/AxisWeightEstimates.lean`
  - `parameterPrimitiveJet`
  - `parameterPrimitiveJet_bound`
- `NavierStokes/AxisOperators.lean`
  - `parameterPrimitiveData`
  - `parameterPrimitive`
  - `norm_parameterPrimitive_le`
- `NavierStokes/AxisEvaluationAlgebra.lean`
  - `parameterPrimitive_eq_deriv_eta`
- `NavierStokes/AxisContraction.lean`
  - `AxisOperators.parameterPrimitive` field and downstream use in the
    coefficient-space natural remainder.

The pinned coefficient identity is

- output radial row `0`: `0`;
- output row `n >= 1`, parameter order `m`:
  `J_out[n,m](eta) = J_in[n-1,m+1](eta) / n`.

Equivalently, this is the coefficient map of the zero-at-axis radial primitive
applied to the parameter derivative,
`integral_0^Y partial_eta F(s,eta) ds`.
The `m+1` shift is part of the official operator; it is not inferred from
numerical differentiation.

## Local artifacts

- `src/openai_ns_reconstruction/axis_coefficient_parameter_primitive.py`
- `tests/test_axis_coefficient_parameter_primitive.py`
- `references/provenance_manifest_addendum_axis_coefficient_parameter_primitive.json`

The production implementation consumes `AxisCoefficientJetState` and obtains
the shifted derivative only from the state's already landed analytic jet
provider. It accepts no caller-supplied derivative table, integration constant,
replacement epsilon, or radial scale.

## Independent checks

The regression suite deliberately does not validate the production row map only
against itself. On the actual schedule-derived angular reference state it:

1. forms a finite radial polynomial from zeroth coefficients;
2. differentiates that polynomial in `eta` by an independent centered finite
   difference;
3. integrates the resulting function in the physical radial variable using
   SciPy adaptive quadrature; and
4. compares that value with the polynomial reconstructed from the production
   parameter-primitive coefficients.

A second finite-difference check targets the next parameter jet, and the actual
axial reference checks the exact degree-one to degree-two relation with the
required `m+1` offset. Index and pinned-window guards remain fail-closed.

## Truth boundary

Status: `formal-structure`.

`paper_exact_velocity_available=false` remains mandatory. This increment does
not certify the global all-index weighted `AxisSpace` norm, instantiate all
remaining linear/composite operators needed by `coefficientOperators`, evaluate
`naturalRemainder(x0)`, produce the first genuine Picard iterate, or materialize
the fixed-point `phi/u`, average/pressure, and `NaturalProfileAssembly`.
Support, moment, matching, and cone closure for the final leading profile also
remain downstream obligations.
