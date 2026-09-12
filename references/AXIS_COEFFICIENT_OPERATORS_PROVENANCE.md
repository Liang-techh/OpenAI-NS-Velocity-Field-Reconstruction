# Axis coefficientOperators provenance

Truth status: **formal-structure**. `paper_exact_velocity_available=false` remains mandatory.

## Source pin

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Primary definition: `NavierStokes/AxisContraction.lean`, `coefficientOperators`.

At the pinned source, the record is assembled as:

- `product := AxisOperators.product I hε`
- `average := AxisOperators.average I hε`
- `primitive := AxisOperators.primitive I hε`
- `parameterPrimitive := AxisOperators.parameterPrimitive I hε`
- `mulY := AxisOperators.mulY I hε`
- `j1/j2 := AxisOperators.regularInverse ... 1/2`
- `param1/param2 := AxisOperators.inverseParamProduct ... 1/2`
- `dot1/dot2 := AxisOperators.inverseDotProduct ... 1/2`
- `mixed1/mixed2 := AxisOperators.inverseMixed ... 1/2`.

## Executable assembly landed here

`src/openai_ns_reconstruction/axis_coefficient_operators.py` groups the already-landed exact jet maps into one `AxisCoefficientOperators` object. The public factory requires `ActualScheduleReferenceAxisState`, inherits its theorem-selected `epsilon`, and rejects operands at any different coefficient scale. It accepts no caller-supplied operator table, radial divisor, integration constant, derivative table, replacement coefficients, or unrelated epsilon.

The record exposes the exact pinned Lean field names `product`, `average`, `primitive`, `parameterPrimitive`, `mulY`, `j1`, `j2`, `param1`, `param2`, `dot1`, `dot2`, `mixed1`, and `mixed2`. Snake-case aliases exist only for Python ergonomics and delegate to the exact fields.

## Verification

`tests/test_axis_coefficient_operators.py` checks all thirteen record fields against their independently landed primitive constructors on the actual SchedulePressure-derived angular/axial reference jets, checks the exact Lean field layout, checks common-epsilon rejection, and confirms every returned state remains fail-closed for global `AxisSpace` membership.

This regression is an assembly test, not a second implementation of the coefficient formulas; each underlying primitive already has its own independent quadrature / finite-difference / physical-identity regression.

## Remaining boundary

This closes the representation-level `coefficientOperators` assembly seam only. It does **not** materialize the global all-index weighted `AxisSpace` norm certificate or continuous-linear/bilinear operator objects from Lean. The next direct Issue #1 blocker is an executable pinned `naturalRemainder` on this record and the actual coefficient data, followed by `naturalRemainder(x0)`, the first genuine Picard iterate, the fixed-point `phi/u`, derived average/pressure, `NaturalProfileAssembly`, and final support/moment/matching/cone checks.

No result in this increment is promoted to paper-exact.
