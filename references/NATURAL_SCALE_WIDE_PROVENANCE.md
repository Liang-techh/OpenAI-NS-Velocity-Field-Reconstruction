# Wide natural-profile Lambda/C selection provenance

Status: **formal-structure**. This increment changes only the scalar representation of the already-landed theorem-side Stage-1 scale selection. It does not construct the coefficient-space fixed point and does not make any leading profile paper-exact.

## Pinned source

Official Lean repository: `openai/NavierStokesAndEuler`

Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant modules:

- `NavierStokes/AxisContraction.lean`
- `NavierStokes/AxisEvaluation.lean`
- `NavierStokes/AxisReference.lean`
- `NavierStokes/NaturalAxisCoefficients.lean`
- `NavierStokes/NaturalProfile.lean`

## Theorem algebra carried in wide form

The existing `natural_scale_selection.py` records the pinned deterministic choices

```text
contractionThreshold = 1 + remainderBound + remainderLip
stabilityScale       = 1 + (14000/9) remainderBound
Lambda               = max(contractionThreshold, stabilityScale)
C                    = exp(Lambda * realPartSup(axisPhase, compactSet)).
```

PR #67 landed `axis_remainder_wide_bounds.py`, which evaluates the same positive `AxisContraction.Controlled` / `controlledRemainder` majorant algebra in upward-rounded 96-digit `Decimal` arithmetic. On the current theorem-admissible SchedulePressure regression (`P=2, m=1, lambda=0.05, wait=30, h=0.01, j=0.05`) at least one conservative remainder majorant is larger than `sys.float_info.max`, so the prior binary64 scale selector cannot consume the certified upper ledger.

`natural_scale_selection_wide.py` now carries exactly the same positive threshold algebra without converting those remainder bounds back to binary64. The rational factor `14000/9` is evaluated with upward rounding. The normalization constant is represented symbolically as

```text
SymbolicExponentialThreshold(exponent_upper)
exponent_upper = Lambda * phase_real_part_sup.
```

Thus `C` means the positive real threshold `exp(exponent_upper)`; the implementation does not attempt to materialize that possibly astronomical exponential in a finite floating-point format. Comparisons of candidate symbolic `C` values are made at the exponent level using monotonicity of `exp`.

`stage1_scale_chain_wide.py` starts from the existing actual-schedule diagnostic, reuses its SchedulePressure analytic neighborhood, componentwise Cauchy ledger, and representable chi-specific AxisResolvent majorant, then evaluates the wide remainder certificate and feeds it directly into the wide Lambda/C selector. It accepts no caller-supplied `rho`, field bounds, resolvent norm, remainder constants, `Lambda`, or `C` on this actual-schedule path.

## Independent regression boundary

Tests cross-check the wide selector against the existing binary64 selector where all quantities are representable, verify the closed rational threshold algebra independently, verify that a symbolic exponential threshold survives values far beyond direct floating-point exponentiation, and run the landed actual SchedulePressure regression through the complete wide remainder -> Lambda/C selection path.

These tests validate the representation bridge only. They do not infer sharper theorem constants from samples and do not solve the coefficient-space fixed-point equations.

## What this increment closes

The former representation blocker after PR #67 is removed: actual SchedulePressure-derived wide `remainderBound/remainderLip` values can now reach a finite wide `Lambda` and a well-defined symbolic theorem-side `C = exp(Lambda * phase_sup)` without narrowing to binary64.

## What remains open

Stage 1 remains **formal-structure** and `paper_exact_velocity_available` must remain `false`. In particular, this increment does not:

- materialize the coefficient-space fixed-point fields `phi/u`;
- derive the corresponding `average/pressure` fields;
- instantiate `NaturalProfileAssembly` with those fields;
- prove support, moments, matching, or cone conditions for the resulting profile;
- turn the conservative Python bounds into interval-arithmetic or Lean proof objects.

The next theorem-faithful blocker is no longer the scalar Lambda/C representation. It is to instantiate the actual coefficient-space fixed-point map with the landed schedule/operator/remainder data (or tighten the conservative bounds if direct numerical iteration requires a more practical representation), then derive `average/pressure` and feed those genuine fields into `NaturalProfileAssembly`.
