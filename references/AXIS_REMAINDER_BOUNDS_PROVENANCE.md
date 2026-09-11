# Axis remainder bound certificate provenance

Status: **formal-structure**. This file does not make Stage 1 paper-exact and does not materialize the missing coefficient-space fixed point.

## Pinned source

- Official Lean repository: `openai/NavierStokesAndEuler`
- Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisContraction.lean`
  - `Controlled`
  - `controlledRemainder`
  - `remainderBound`
  - `remainderLip`
  - `contractionThreshold`
  - `coefficient_fixedPoint`
- `NavierStokes/AxisOperators.lean`
  - `norm_product_le`: `64`
  - `norm_average_le`: `1`
  - `norm_primitive_le`: `80`
  - `norm_regularInverse_le`: `80`
  - `norm_parameterPrimitive_le`: `80 / epsilon`
  - `norm_mulY_le`: `80`
  - `norm_inverseParamProduct_le`: `5120 / epsilon`
  - `norm_inverseDotProduct_le`: `5120`
  - `norm_inverseMixed_le`: `5120 / epsilon`
- `NavierStokes/AxisResolvent.lean`
  - `naturalResolvent_norm_le`: the resolvent norm is bounded by the factorial-majorant series and is therefore kept as an independently certified upstream input in this increment.

## Executable mapping

`src/openai_ns_reconstruction/axis_remainder_bounds.py` mirrors the two scalar fields propagated by Lean's `Controlled` constructors. Given certified upper bounds on the fixed coefficient norms and operator norms, it computes conservative upper bounds for the complete `controlledRemainder` bound and Lipschitz fields. The module also reconstructs a conservative radius from the actual reference-pair formula

`(S one, -(1/2) * j1(product inverseL zStar))`

using Mathlib's product max norm before running the remainder bookkeeping.

`NaturalOperatorNormBounds.pinned_axis_operators(...)` encodes only operator inequalities explicitly proved in the pinned Lean source. It deliberately requires `resolvent_norm_upper` from upstream instead of replacing the factorial series by an arbitrary fitted constant.

The derived coefficient norm bounds use triangle inequality plus the certified product-operator norm for exactly the six coefficient combinations appearing in `naturalRemainder`:

- angular linear coefficient,
- angular quadratic coefficient,
- average coefficient,
- angular slow coefficient,
- axial linear coefficient,
- axial quadratic coefficient.

## What this closes

The previous Stage-1 blocker described `remainderBound/remainderLip` as opaque caller-supplied scalars. This increment replaces that opaque interface with a deterministic theorem-shaped propagation from a finite norm ledger. Once the actual `AxisData` coefficient norms, normalized-amplitude norm `M`, and the factorial-series resolvent norm are independently certified, the existing `natural_scale_selection.py` can consume these conservative bounds to choose a safe `Lambda` without sampling or fitting the nonlinear remainder.

## What remains open

This is **not** a certificate for the actual coefficient family yet. The following inputs still need to be materialized from the paper/Lean construction:

1. the true coefficient-space `AxisData` elements and certified norms (`one`, `eta`, `d`, `inverseL`, `uStar`, `uStarEta`, `wStar`, `hStar`, normalized gradient, `zStar`);
2. the actual normalized angular input bound `M`;
3. a certified numerical upper bound for the `naturalResolvent` factorial-majorant series from the actual `chi` norm;
4. the complex compact-set `realPartSup` used by the normalization threshold;
5. the actual fixed-point pair `phi/u`, followed by `average/pressure` and `NaturalProfileAssembly`;
6. support, moment, matching, and cone-condition verification for the materialized profile.

Therefore `paper_exact_velocity_available` must remain `false`, and Stage 1 remains `formal-structure`.

## Test boundary

`tests/test_axis_remainder_bounds.py` checks the pinned operator constants, independent algebra for the six derived coefficient bounds, the reference-pair max-norm radius, an exact scalar `Controlled` wiring case, monotonicity in the normalized amplitude bound, and fail-closed invalid inputs. The simplified scalar wiring case is only a unit test of bookkeeping and is explicitly not a surrogate Navier--Stokes profile.
