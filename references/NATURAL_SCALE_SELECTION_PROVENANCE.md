# Natural-profile Lambda/C selection provenance

Status: **formal-structure**. This file records a theorem-faithful parameter-selection interface. It does **not** claim that the coefficient-space constants feeding the interface have yet been materialized from the actual Stage 1 analytic coefficient family.

## Pinned source

Official Lean repository: `openai/NavierStokesAndEuler`

Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant modules:

- `NavierStokes/AxisContraction.lean`
- `NavierStokes/AxisEvaluation.lean`
- `NavierStokes/AxisReference.lean`
- `NavierStokes/NaturalAxisCoefficients.lean`
- `NavierStokes/NaturalProfile.lean`

## Exact selection chain implemented

`AxisContraction.contractionThreshold` is

```text
1 + remainderBound + remainderLip
```

where both controlled-remainder constants are evaluated on the radius `||referencePair|| + 1` and the uniform angular-amplitude norm bound `M`. `uniform_natural_fixedPoint` uses this quantity directly as its fixed-point scale threshold.

`AxisReference.exists_positive_scaled_profiles` then replaces that threshold by the maximum with

```text
stabilityScale epsilon K
  = 1 + 500 * (jetBound epsilon 5 0 0 + jetBound epsilon 5 1 0) * K,
```

where `K = errorConstant`, and `NaturalProfile.profileErrorConstant` is the same evaluated `remainderBound` occurring in the coefficient-space error estimate.

The two required `AxisEvaluation.jetBound` values do not need to remain caller inputs. From the pinned `majorant` definition, when `R=5` and `m=0`, the epsilon-dependent prefactor is exactly one. Hence

```text
jetBound epsilon 5 0 0
  = sum_n (1/4)^n
  = 4/3,

jetBound epsilon 5 1 0
  = sum_n (n+1) (1/4)^n
  = 16/9.
```

Therefore the landed stability threshold closes to

```text
stabilityScale epsilon B = 1 + (14000/9) B.
```

Consequently, after the actual coefficient-space remainder constants have been evaluated, `natural_scale_selection.py` selects

```text
Lambda = max(
    1 + B + L,
    1 + (14000/9) * B,
)
```

with `B = remainderBound` and `L = remainderLip`.

Finally, `NaturalAxisCoefficients.AnalyticInputs.normalizationThreshold` is literally

```text
exp(Lambda * realPartSup(axisPhase, compactSet)),
```

and `NaturalProfile.exists_natural_profiles` chooses `C` at equality with that threshold. The runtime interface therefore sets

```text
C = exp(Lambda * phase_real_part_sup).
```

This is not a tunable heuristic and no sampled phase maximum is substituted for `realPartSup`.

## What this increment closes

The repository now has a single executable record of the exact theorem-side algebra that chooses `Lambda` and `C` once the actual coefficient-family remainder constants and compact complex phase supremum are available. The two `R=5` jet sums are analytically closed rather than caller supplied. The module also has a fail-closed admissibility checker for larger candidate `Lambda,C` pairs, recomputing the normalization threshold at the candidate `Lambda`.

## What remains open

The following are still required before Stage 1 can be paper-exact:

- construct the actual common analytic coefficient family for the schedule-derived pressure datum in executable coefficient space;
- evaluate/certify its `remainderBound` and `remainderLip` rather than supplying them to this adapter;
- evaluate/certify `realPartSup(axisPhase, compactSet)` on the actual compact complex neighborhood rather than supplying it;
- use those certified constants to materialize the coefficient-space fixed point `phi/u`, derive `average/pressure`, and feed the resulting fields into `NaturalProfileAssembly`;
- verify the final support, moment, matching, and cone conditions.

Accordingly, this increment must not change `paper_exact_velocity_available` or Stage 1 from `formal-structure`.

## Numerical boundary

The theorem works over real numbers. The Python adapter stores finite binary64 values and refuses non-finite/negative theorem constants. If `exp(Lambda * phase_real_part_sup)` exceeds binary64 range it fails closed instead of returning an infinite normalization and pretending the construction remains executable. Tests verify the closed jet sums, selection identities and inequality gate; they do not certify that caller-supplied `B`, `L`, or phase supremum are the actual coefficient-space quantities.
