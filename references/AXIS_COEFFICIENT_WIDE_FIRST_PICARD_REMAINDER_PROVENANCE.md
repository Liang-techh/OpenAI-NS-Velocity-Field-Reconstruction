# Complete first-Picard natural remainder provenance

Status: **formal-structure only**.

Artifact:

- `src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_remainder.py`

This increment assembles the complete coefficientwise pinned
`naturalRemainder(x1)` from the actual theorem-selected first Picard state.
The angular half is

```text
naturalResolvent(inverseL * (lin1 + quad1 - slow1/Lambda))
```

and the axial half is

```text
inverseL * (lin2 - slow2/Lambda + pressure).
```

The raw constituent modules are evaluated on the exact same typed `x1`.
Scale powers are retained as normalized Decimal numerators: ordinary channels
`Lambda^0` through `Lambda^-5`, pressure-linear channels `a^2 Lambda^0`
through `a^2 Lambda^-4`, and the axial pressure-square channel
`a^4 Lambda^-3`.  The angular pressure-square channel is structurally zero.
The pressure state supplies full eta derivatives of the normalized pressure,
including amplitude derivatives; the remainder assembly does not differentiate
that value again.

The angular natural resolvent is evaluated coefficientwise by its finite radial
filtration.  With `chi` radially constant, its recurrence is

```text
R[n,m] = A[n,m]
       - sum(k=0..m) binom(m,k) chi[0,k] R[n-1,m-k]
         / (2*n*(n+1)),       n >= 1,
R[0,m] = A[0,m].
```

The outer axial expression stops after `inverseL`; it does not receive this
angular resolvent.  Per-request memoization caches source and resolved rows so
the finite recurrence does not repeat constituent evaluations.

The same module exposes a coefficientwise second Picard view,
`referencePair + naturalRemainder(x1)/(2*Lambda)`.  Its reference channel is
read directly from the genuine actual reference state.  Remainder numerators
are shifted one scale power and divided by two, yielding ordinary channels
through `Lambda^-6`, pressure-linear channels through `a^2 Lambda^-5` with a
structural zero at power zero, and normalized pressure-square `a^4 Lambda^-4`.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main files:
  - `NavierStokes/AxisContraction.lean`
  - `NavierStokes/AxisOperators.lean`
  - `NavierStokes/AxisWeightEstimates.lean`
  - `NavierStokes/NaturalAxisCoefficients.lean`
- Relevant symbols:
  - `NavierStokes.AxisContraction.naturalRemainder`
  - `NavierStokes.AxisResolvent.naturalResolvent`
  - `NavierStokes.AxisOperators.regularInverse`
  - `NavierStokes.AxisOperators.inverseDotProduct`
  - `NavierStokes.AxisOperators.inverseMixed`
  - `NavierStokes.AxisOperators.inverseParamProduct`
  - `NavierStokes.AxisOperators.average`
  - `NavierStokes.NaturalAxisCoefficients.CoefficientFamily.axisData`

The coefficient definitions and natural-remainder formula are at
`AxisContraction.lean:327-365`.  Operator bindings are at
`AxisOperators.lean:86-95`, `AxisOperators.lean:501-518`, and
`AxisOperators.lean:681-717`; the radial divisor is pinned at
`AxisWeightEstimates.lean:388`; and the parameter window is pinned at
`NaturalAxisCoefficients.lean:27`.

## Domain and truth boundary

The production factory accepts only a genuine
`ActualScheduleWideFirstPicardState`.  It rebuilds fixed `AxisData` and `chi`
from `x1.reference`, checks the actual epsilon, certified schedule data, and
positive theorem-selected `Lambda`, and requires every constituent branch to
refer to that identical `x1`.  Coefficient indices satisfy `n,m >= 0` and
`eta` lies in the pinned window `[-11/10, 11/10]`.

`natural_remainder_x1_materialized` and the coefficientwise second-Picard
marker are true after all constituent APIs pass their checks.  The result
keeps `paper_exact`, `global_axis_norm_certified`, `fixed_point_materialized`,
and `fixed_point_convergence_certified` false.  This boundary does not claim
global weighted `AxisSpace` membership, infinite-series convergence, a final
fixed point, downstream derived fields, or a paper-exact velocity profile.

The Decimal arithmetic uses a local precision of 96 digits.  It is a numerical
representation of the pinned finite coefficient formulas, not exact real
arithmetic, an interval enclosure, or an independent theorem certificate.
