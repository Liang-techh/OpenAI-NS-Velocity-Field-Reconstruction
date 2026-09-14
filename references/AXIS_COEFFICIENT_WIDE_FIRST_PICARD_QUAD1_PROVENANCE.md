# Genuine first-Picard `quad1(x1)` provenance

## Scope

This increment materializes exactly the angular quadratic branch

`quad1(x1) = -j2(product(primitive(phi1), phi1))`

on the landed genuine `ActualScheduleWideFirstPicardState`.  The first Picard
angular state is not projected to binary64; its `Lambda^0`, `Lambda^-1`, and
`Lambda^-2` coefficient families are propagated independently through the
pinned primitive/product/j2 rules, producing `Lambda^0` through `Lambda^-4`.

It does **not** claim that `naturalRemainder(x1)` is complete, does not produce
`x2`, and does not certify a global `AxisCoefficientSpace` norm, contraction
closure, a fixed point, `NaturalProfileAssembly`, or paper-exact velocity.

## Formal source

Pinned source revision:

- repository: `openai/NavierStokesAndEuler`
- commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/AxisContraction.lean`
  - `NavierStokes.AxisContraction.naturalRemainder`
  - literal branch `quad1 := - O.j2 (O.product (O.primitive phi) phi)`
- `NavierStokes/AxisOperators.lean`
  - `NavierStokes.AxisOperators.primitive`
  - `NavierStokes.AxisOperators.product`
  - `NavierStokes.AxisOperators.regularInverse`
  - `primitiveScale`

The replayed row identities are:

- `primitive(phi)[0,m] = 0`;
- `primitive(phi)[n,m] = phi[n-1,m]/n` for `n>0`;
- product uses the full radial convolution and eta-Leibniz sum;
- `j2(f)[0,m] = 0`;
- `j2(f)[n,m] = f[n-1,m]/(n(n+1))` for `n>0`.

The sole minus sign is the literal one in the `quad1` definition.

## Actual-data dependency

The only field input is the same theorem-selected
`ActualScheduleWideFirstPicardState` already produced from the actual
SchedulePressure reference chain.  Its angular coefficient jet is

`phi1 = r + b/Lambda + c/Lambda^2`.

No coefficient table, cutoff, Lambda, amplitude, surrogate field, fitted
sample, or omitted-row value is accepted from the caller.  The only explicit
zero rows are the pinned row-zero semantics above.

## Regression obligations

`tests/test_axis_coefficient_wide_first_picard_quad1.py` checks:

1. source-revision-bound genuine `x1` typing and fail-closed truth flags;
2. the literal `j2` row-zero identity;
3. the first nonzero radial row identity
   `quad1[2,0] = -phi1^2[0,0]/6` coefficient-family by coefficient-family;
4. a general `(n,m)` jet against an independently replayed
   primitive/product/j2 finite sum using exact `Decimal` arithmetic;
5. preservation of every inverse-Lambda family without an O(1) scale collapse;
6. rejection of surrogate states and invalid indices/window coordinates; and
7. machine-readable provenance with `naturalRemainder(x1)`, `x2`, fixed-point,
   and paper-exact claims still false.

## Remaining boundary

The genuine x1 source, pressure, complete four-constituent `slow2(x1)`, and now
`quad1(x1)` are materialized on compatible actual-schedule data.  The shortest
remaining closure path is to materialize `lin1(x1)`, `slow1(x1)`, and `lin2(x1)`
in compatible mixed-scale typed form, then perform the complete angular/axial
`naturalRemainder(x1)` recombination.  Only after that recombination may the
outer Picard map construct a genuine `x2`.  Iteration, contraction/tail closure,
global coefficient-space membership, fixed-point `phi/u`, derived
average/pressure fields, `NaturalProfileAssembly`, and paper-exact velocity
remain unresolved.
