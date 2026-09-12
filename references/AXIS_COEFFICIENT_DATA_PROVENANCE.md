# Axis coefficient data provenance

This increment materializes the fixed coefficient fields required by the pinned
`AxisContraction.AxisData` / `NaturalAxisCoefficients.CoefficientFamily.axisData`
record on the same actual SchedulePressure-derived coefficient scale already used
by the landed reference pair and `coefficientOperators`.

Source pin: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The official `CoefficientFamily.axisData` record supplies the scalars
`A = NaturalAxisData.A h`, `D = NaturalAxisData.D h`, `h`, and the radially
constant fields `one`, `eta`, `d`, `inverseL`, `uStar`, `uStarEta`, `wStar`,
`hStar`, `normalizedGradient`, and `zStar`.  The coefficient-family theorem
states that these fields have only radial degree zero.  The Python realization
therefore sets every radial row `n>0` exactly to zero and computes arbitrary
finite eta derivatives of row zero by normalized truncated-Taylor algebra.

For `zStar`, the implementation uses the actual landed
`SchedulePressure.axisPressure` Taylor chain via
`axis_pressure_normalized_taylor`; it does not substitute the earlier ideal
prefix pressure witness, fit samples, or accept caller-supplied pressure data.
The factory accepts only `ActualScheduleReferenceAxisState`, so `h`, `j`,
`sigma`, and `epsilon` are inherited from the same theorem-selected schedule and
analytic-neighborhood chain as the reference pair.

Independent regression checks row-zero values against the direct real-field
formulas in `natural_axis.py`, including a direct `zStar` evaluation from the
actual `axis_pressure` and `axis_pressure_derivative` evaluators.  First eta
jets for representative rational/polynomial/pressure-dependent fields are
cross-checked against centered finite differences of those independent point
functions.  The assembled `coefficientOperators` record is also required to
accept every fixed field without an epsilon mismatch.

## Truth boundary

This is a representation-level materialization only.  It does **not** certify
the global all-index weighted `AxisSpace` norm of the fixed fields, instantiate
the Lean continuous-linear/bilinear norm objects, execute the pinned
`naturalRemainder`, evaluate `naturalRemainder(x0)`, produce a genuine Picard
iterate or fixed point, derive the final average/pressure fields, connect those
fields to `NaturalProfileAssembly`, or close support/moment/matching/cone
conditions.  Consequently Stage 1 remains `formal-structure` and
`paper_exact_velocity_available=false`.
