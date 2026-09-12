# Stage-1 mixed-scale angular natural remainder provenance

Status: **formal-structure**. This increment does not claim a paper-exact leading profile or a materialized fixed point.

## Pinned source

Mapped to `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, principally:

- `NavierStokes/AxisContraction.lean`
  - `naturalRemainder`
  - the angular `lin1`, `quad1`, `slow1` branch
- `NavierStokes/AxisOperators.lean`
- `NavierStokes/AxisResolvent.lean`
- `NavierStokes/NaturalAxisCoefficients.lean`

At the contraction centre `x0 = referencePair`, with theorem-selected
`t = 1/Lambda`, the pinned angular remainder is

`naturalResolvent(inverseL * (lin1 + quad1 - slow1 / Lambda))`.

## Landed capability

`src/openai_ns_reconstruction/axis_coefficient_wide_angular_remainder.py`
binds this formula to the actual `SchedulePressure` reference state, the landed
13-field coefficient operator record, the actual `AxisData`, and the
coefficientwise-exact natural resolvent.

The implementation uses linearity to retain the theorem-selected scale without
narrowing `1/Lambda` to binary64:

1. `naturalResolvent(inverseL * (lin1 + quad1))` is evaluated as the ordinary base state;
2. `naturalResolvent(inverseL * slow1)` is evaluated as an ordinary numerator state;
3. each coefficient jet is represented as `ordinary_base - ordinary_slow/Lambda`, with the actual theorem-selected Decimal `Lambda` kept explicitly.

This is exactly the same pinned expression because multiplication by the fixed
`inverseL` field and `naturalResolvent` are linear. The natural resolvent itself
uses the already-landed radial filtration identity: row `n` is the exact finite
alternating sum through `Q^n`, not a tolerance-based truncation.

There is no caller override for sigma, epsilon, `Lambda`, `C`, coefficient
tables, resolvent cutoff, or derivative tables. `Lambda` is obtained from the
same actual SchedulePressure wide scale chain used by the landed theorem
amplitude.

## Independent verification

The regression suite cross-checks both ordinary numerator states against the
previously landed full executable `naturalRemainder` composition:

- regression-only `t=0, a=0` isolates
  `naturalResolvent(inverseL * (lin1 + quad1))`;
- with the same `a=0`, the difference between the angular components at
  regression-only `t=0` and `t=1` isolates
  `naturalResolvent(inverseL * slow1)`.

Those values are test oracles only and are never promoted to manuscript scale
choices. The tests also require the production `Lambda` to agree with the
actual theorem-selected amplitude scale, require a nonzero slow numerator to
remain a nonzero Decimal after division by `Lambda`, and check fail-closed
index/window handling.

## Explicit non-claims / next blocker

This closes the **angular** mixed-scale branch only. The separately landed axial
branch still carries the genuine signed-log `a^2` pressure scale. The next
smallest blocker is to assemble these two typed decompositions as the complete
mixed-scale pair `naturalRemainder(x0)` and then apply the outer fixed-point
factor `1/(2*Lambda)` without scale collapse to materialize the first genuine
Picard iterate `x1`.

No fixed-point convergence/certificate, final `phi/u`, radial average/pressure
reconstruction, `NaturalProfileAssembly`, global all-index weighted `AxisSpace`
certificate, or support/moment/matching/cone closure is claimed.

`paper_exact_velocity_available=false` remains mandatory.
