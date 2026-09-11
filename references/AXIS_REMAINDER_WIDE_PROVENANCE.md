# Axis remainder wide-arithmetic provenance

Status: **formal-structure / representation-layer certificate only**.

This increment does not materialize the coefficient-space fixed point and does
not change `paper_exact_velocity_available=false`.

## Pinned source

Official source commit:
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The scalar algebra being mirrored is the same one already mapped in
`AXIS_REMAINDER_BOUNDS_PROVENANCE.md` to the pinned natural-axis contraction
construction, in particular `AxisContraction.Controlled` and
`controlledRemainder`, together with the explicit operator norm inequalities
from `AxisOperators.lean`.

The previous executable adapter
`src/openai_ns_reconstruction/axis_remainder_bounds.py` stores every propagated
bound as binary64.  For the actual landed SchedulePressure analytic ledger,
the chi-specific AxisResolvent majorant is now representable after PR #62, but
positive products in the subsequent conservative `remainderBound` /
`remainderLip` bookkeeping exceed binary64 range.

## What this increment adds

`src/openai_ns_reconstruction/axis_remainder_wide_bounds.py` evaluates the same
positive triangle/product/Lipschitz majorant algebra with `decimal.Decimal`.
Every finite binary64 input is imported exactly with `Decimal.from_float`, and
positive additions/multiplications performed by this module are rounded toward
`+infinity` at 96 decimal digits.  Consequently this layer does not respond to
overflow by clipping, rescaling, fitting, or substituting smaller sampled
values.

The adapter includes wide versions of:

- the derived natural-axis coefficient norm ledger;
- the reference-pair max-norm bound and radius `||referencePair||+1`;
- the `Controlled` constant/coordinate/add/sub/linear/bilinear rules;
- the full natural-axis remainder expression, including pressure terms; and
- the final conservative `remainderBound` and `remainderLip` values.

Tests cross-check an exact simple `Controlled` case, compare against the
binary64 implementation when all intermediates are representable, and then run
the already-landed actual SchedulePressure regression datum
`P=2, m=1, lambda=0.05, wait=30, h=0.01, j=0.05`.  On that actual datum the
wide remainder ledger stays finite in Decimal even though at least one final
remainder majorant is larger than `sys.float_info.max`.

This identifies the PR #62 blocker as a **binary64 representation boundary in
the conservative upper ledger**, not a failure to propagate the theorem-shaped
remainder estimate and not a lower bound on the true remainder.

## Truth boundary and next step

The upstream analytic/operator inputs remain executable real-arithmetic
certificates rather than interval or Lean proof objects.  Wide Decimal
arithmetic only preserves the direction of the supplied positive majorant
algebra; it does not upgrade those upstream inputs to formal proof.

This increment still does **not**:

- produce Banach-space elements `phi` or `u`;
- derive `average` or `pressure` fixed-point fields;
- instantiate `NaturalProfileAssembly`;
- prove support, moment, matching, or cone conditions; or
- make the enormous wide remainder constants numerically useful by itself.

The next smallest theorem-faithful step is to carry these wide
`remainderBound/remainderLip` values through the pinned `Lambda/C` scale
selection without converting them back to binary64, and only then determine
whether a tighter theorem-side majorant is required before an actual
coefficient-space fixed-point iteration can be instantiated.
