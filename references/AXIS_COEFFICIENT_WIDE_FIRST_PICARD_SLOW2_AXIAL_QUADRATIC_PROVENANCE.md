# Provenance: theorem-scale first-Picard `slow2` axial quadratic branch

## Truth status

- Layer: `stage-1-leading-profile`
- Status: `formal-structure`
- `paper_exact_velocity_available = false`
- `full_reconstruction = false`

This artifact closes one narrow coefficient-space seam only. It does not claim complete `slow2(x1)`, complete `naturalRemainder(x1)`, `x2`, fixed-point convergence/finality, global weighted `AxisSpace` membership, `NaturalProfileAssembly`, or paper-exact leading velocity.

## Pinned source

Official Lean repository and revision:

- `openai/NavierStokesAndEuler`
- commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant definitions:

- `NavierStokes/AxisContraction.lean`: `naturalRemainder`, whose first `slow2` term is `j1 (product axialQuadraticCoefficient (product u u))`.
- `NavierStokes/AxisContraction.lean`: `axialQuadraticCoefficient = 2 * A * eta`.
- `NavierStokes/AxisWeightEstimates.lean` / `NavierStokes/AxisOperators.lean`: pinned coefficient product / eta-Leibniz action.
- `NavierStokes/AxisOperators.lean`: zero-datum `regularInverse`; for `j1`, output row `n>=1` reads source row `n-1` divided by `radialDivisor(1,n-1)=n^2`, and row zero vanishes.

## Actual-data chain

Input is only the already-landed genuine `ActualScheduleWideFirstPicardState`.

The state therefore inherits the actual SchedulePressure datum, selected `j`, certified `sigma`, coefficient-space `epsilon`, wide theorem-selected `Lambda`, and genuine signed-log amplitude/pressure chain. The scalar `A` is recovered from `actual_schedule_axis_coefficient_data(x1.reference)`; callers cannot replace it.

The input `u1*u1` is the landed theorem-scale split product from `axis_coefficient_wide_first_picard_axial_square.py`, retaining separately:

- ordinary `Lambda^0 ... Lambda^-4` numerators;
- pressure-linear `a^2 Lambda^-1 ... Lambda^-3` numerators;
- pressure-square `a^4 Lambda^-2` numerator.

No toy amplitude, fitted coefficient/derivative table, caller-selected `Lambda/C`, pressure replacement, or arbitrary radial/product cutoff is used.

## Executed identity

Let `F = u1*u1` and `q(eta)=2*A*eta`. Because `q` is radial-degree zero, `q'=2*A`, and higher eta derivatives vanish,

`d_eta^m(q F) = q d_eta^m F + m (2*A) d_eta^(m-1) F`.

This identity is applied independently to every retained theorem-scale numerator. The pinned `j1` is then applied exactly:

- output row `0` is zero;
- output row `n>=1` is the transformed source row `n-1` divided by `n^2`.

The common signed-log `a^2` and `a^4` scales are not converted to binary64. Only normalized Decimal numerators are transformed, so mathematically nonzero theorem-scale terms cannot silently underflow to zero.

## Validation

Regression coverage verifies:

1. The constructor accepts only the genuine theorem-scale first Picard state and keeps all final/fixed-point flags false.
2. The ordinary `Lambda^0` component agrees with an independent composition of the previously landed pinned operators: `u0*u0 -> eta*(u0*u0) -> j1 -> 2*A`.
3. Every ordinary, `a^2`, and `a^4` numerator agrees exactly with a separate literal implementation of the eta product rule followed by division by `n^2`.
4. The `j1` output row zero is exactly zero while still routing through the actual upstream state guards.
5. Nonzero pressure-linear terms retain common log scale `2*log(a)` and the pressure-square term retains `4*log(a)`, with explicit inverse-`Lambda` powers only in the signed-log factor.
6. Current nonzero theorem-scale pressure pieces fail closed on binary64 underflow rather than becoming zero.
7. Invalid indices, out-of-window eta, and surrogate/non-first-Picard state inputs fail closed.

## Remaining boundary

The other `slow2(x1)` terms still require mixed-scale propagation, in particular the branches coupling `average(u1)` with `u1` through `dot1`, `mixed1`, and `param1`. The angular `lin/quad/slow` branches and the new source/pressure/resolvent recombination also remain before complete `naturalRemainder(x1)` can be formed. Only after that can a genuine `x2` be materialized and the pinned contraction/convergence argument be attached to obtain the fixed point `phi/u`. Derived average/pressure, `NaturalProfileAssembly`, global all-index weighted `AxisSpace` certification, and support/moment/matching/cone closure remain downstream.
