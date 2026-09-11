# Stage-1 natural-axis Picard bridge provenance

Status: **formal-structure**. This note does not claim a paper-exact leading profile.

## Pinned source

This increment is mapped to `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, principally:

- `NavierStokes/AxisContraction.lean`
  - `exists_fixedPoint_of_controlled`
  - `exists_unique_natural_fixedPoint`
  - `coefficient_fixedPoint`
  - `natural_axis_profiles`
- `NavierStokes/AxisCoefficientSpace.lean`
- `NavierStokes/AxisOperators.lean`
- `NavierStokes/AxisResolvent.lean`

The fixed-point map in the pinned theorem is

`F(x) = x0 + (1 / (2 * Lambda)) • naturalRemainder O d S (1 / Lambda) a x`,

with `x0 = referencePair O d S`. The proof of
`exists_fixedPoint_of_controlled` checks the two scalar gates

`(1 / (2 * Lambda)) * remainderBound <= 1`

and

`(1 / (2 * Lambda)) * remainderLip <= 1/2`.

## Landed capability

`src/openai_ns_reconstruction/axis_fixed_point_picard.py` consumes the already
landed **actual SchedulePressure** wide chain; it does not accept independent
caller choices for `remainderBound`, `remainderLip`, or `Lambda` on that path.
It keeps the quantities in 96-digit `Decimal`, forms an upward-rounded upper
bound for `1/(2*Lambda)`, and then checks conservative upper bounds for both
fixed-point gates. This is the first executable bridge from the actual schedule
`rho / coefficient bounds / resolvent / remainder / Lambda` chain to the
specific contraction hypotheses used by the pinned coefficient-space theorem.

The same module exposes the exact Picard-map *shape*. Scalar multiplication by
`1/(2*Lambda)` is deliberately delegated to the future coefficient-space
backend so an enormous theorem-side `Lambda` never has to be narrowed through
binary64 merely to execute the iteration.

The geometric-tail helper is the standard Banach estimate
`q^n/(1-q) * ||x1-x0||`, with `q` replaced by the conservative certified upper
bound produced above. It is useful only after the genuine coefficient backend
can supply the corresponding coefficient-space norm of the first step.

## Verification

The regression suite checks that the existing theorem-admissible schedule
`P=2, m=1, lam=0.05, wait=30, h=0.01, j=0.05` reaches the certificate and
satisfies both scalar gates despite the conservative quantities exceeding
binary64 range. Separate tests cross-check the geometric-tail scalar against an
independent Decimal computation and verify the generic iteration adapter on a
small Decimal algebra fixture. That fixture tests plumbing only; it is not a
surrogate leading profile.

## Explicit non-claims / next blocker

This increment **does not** implement `AxisCoefficientSpace.AxisSpace`, the
pinned `coefficientOperators`, `AxisResolvent.naturalResolvent`, or the actual
`naturalRemainder` evaluator. It therefore does not yet materialize the
coefficient-space fixed point `x=(phi,u)`, its radial average, or its pressure
field. It also does not instantiate `NaturalProfileAssembly`, verify support or
moments, or close matching/cone conditions.

The next theorem-faithful step is to give Python an executable representation
of the genuine compatible coefficient state and enough of the pinned
`coefficientOperators`/`naturalRemainder` to evaluate one real Picard step from
`referencePair`. Once that exists, this certificate supplies the already-fixed
`Lambda` and the contraction/tail budget without introducing a new arbitrary
parameter.
