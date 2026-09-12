# Wide first-Picard coefficient state provenance

## Scope

This artifact records one narrow Stage-1 increment for Issue #1: the first
Picard iterate of the pinned natural-axis contraction is materialized at the
actual SchedulePressure contraction centre without collapsing theorem-selected
wide scales.

The pinned source is `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, primarily
`NavierStokes/AxisContraction.lean`.  The fixed-point step used there has the
form

`F(x) = x0 + (1 / (2 * Lambda)) * naturalRemainder(x)`

with `x0 = referencePair`.  The landed PR #197 path already supplies the
complete mixed-scale `naturalRemainder(x0)` pair on the actual SchedulePressure
construction.  This increment applies only the outer factor and adds `x0`.

## Representation

The theorem-selected `Lambda` is a finite Decimal value far outside the scale at
which a naive binary64 inverse should be used.  Therefore `x1` is retained as a
sum decomposition rather than rounded into a single number.

For one angular coefficient jet,

`phi1 = phi0 + base/(2*Lambda) + slowNumerator/(2*Lambda^2)`.

For one axial coefficient jet,

`u1 = u0 + base/(2*Lambda) + slowNumerator/(2*Lambda^2)
          + pressure/(2*Lambda)`.

The two ordinary correction numerators remain Decimal.  The pressure term
retains the existing split signed-log representation: its common `a^2` log
scale is unchanged and `log(2*Lambda)` is subtracted from the moderate
`log_factor`.  This is deliberate; subtracting that moderate correction from
the current enormous common amplitude log at fixed Decimal precision could
erase it.

The reference value and correction pieces are also deliberately not summed at
fixed precision.  On the current conservative theorem scale an O(1) reference
plus an O(1/Lambda) nonzero correction would otherwise be reported as exactly
the reference value.

## Actual-data boundary

Production construction accepts only the same `TailData`, schedule `j`, and
optional quadrature resolution already used by the landed actual-schedule
chain.  It does not accept overrides for pressure data, sigma, epsilon, Lambda,
C, amplitude, coefficient tables, derivative tables, or radial/resolvent
cutoffs.  The complete remainder object enforces shared TailData, j, certified
sigma, epsilon, and theorem-selected Lambda before this step is reached.

No toy amplitude or arbitrary parameter is promoted to construction data.
Regression values remain regression-only.

## Validation

`tests/test_axis_coefficient_wide_first_picard.py` checks that:

- the state is anchored to the complete actual-schedule `naturalRemainder(x0)`
  and its exact theorem-selected Decimal Lambda;
- angular and axial ordinary pieces satisfy the algebraic outer
  `1/(2*Lambda)` scaling, including the inherited slow `1/Lambda` term becoming
  a `1/Lambda^2` correction;
- the axial signed-log pressure keeps its sign and common amplitude log while
  acquiring exactly the `-log(2*Lambda)` factor correction;
- exact zero pressure remains zero under outer scaling;
- any nonzero pressure that cannot be projected to binary64 fails closed rather
  than becoming zero; and
- fixed-point convergence/finality and paper-exact status remain false.

## Truth boundary

This is a genuine first Picard iterate of the currently materialized
actual-schedule coefficient construction, but it is **not** the fixed point and
it is **not** a paper-exact velocity reconstruction.  The underlying
`realPhase` value still uses the landed numerical quadrature, and the global
all-index weighted `AxisSpace` membership/norm certificate remains absent.

The next Stage-1 blocker is to extend the mixed-scale arithmetic from `x1` to
successive Picard states and supply a convergence/fixed-point certificate for
`phi/u`.  Only after that should derived average/pressure be materialized,
connected to `NaturalProfileAssembly`, and followed by support, moment,
matching, cone-condition, and global `AxisSpace` closure.
