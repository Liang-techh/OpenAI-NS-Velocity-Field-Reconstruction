# Provenance: theorem-scale first-Picard `slow2` average/u `dot1` branch

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

- `NavierStokes/AxisContraction.lean`: `naturalRemainder`, whose `slow2` contains `dot1 (product averageCoefficient bu) u` after `bu = average u`.
- `NavierStokes/AxisContraction.lean`: `averageCoefficient = 2 * D * eta`.
- `NavierStokes/AxisOperators.lean`: `averageData`, giving radial row `n` the exact factor `1/(n+1)`.
- `NavierStokes/AxisOperators.lean`: `dot1 = inverseDotProduct 1 0`. For output row `n>=1`, `i+j=n-1` contributes the radial factor `j / radialDivisor(1,n-1) = j/n^2`; row zero vanishes.
- `NavierStokes/AxisWeightEstimates.lean` / `AxisOperators.lean`: the pinned eta-Leibniz binomial action for coefficient products and bilinear operators.

## Actual-data chain

Input is only the already-landed genuine `ActualScheduleWideFirstPicardState`. It inherits the actual SchedulePressure datum, selected `j`, certified `sigma`, coefficient-space `epsilon`, theorem-selected wide `Lambda`, and genuine signed-log amplitude/pressure chain.

The first Picard axial field is retained as

`u1 = U0 + U1/Lambda + U2/Lambda^2 + a^2 P/Lambda`.

The average is applied componentwise using the exact row factor `1/(n+1)`. The scalar `D` is recovered only from `actual_schedule_axis_coefficient_data(x1.reference)`; callers cannot replace it. Multiplication by `q(eta)=2*D*eta` uses

`d_eta^m(q F) = q d_eta^m F + m*(2*D) d_eta^(m-1) F`.

The common amplitude is the landed theorem-selected `a(eta)=exp(Lambda*realPhase(eta))/C`. Pressure normalized factors already contain the exact eta-jet action of that amplitude before the common `a^2` factor is split, so the product rule above preserves the actual amplitude derivatives rather than using a constant-amplitude surrogate.

No toy amplitude, fitted coefficient/derivative table, caller-selected `D`, `Lambda/C`, pressure replacement, or arbitrary radial/product cutoff is used.

## Executed `dot1` identity

Let `B=(2*D*eta)*average(u1)`. For output row `n>=1`, the implementation executes exactly

`dot1(B,u1)[n,m] = sum_{i+j=n-1} sum_{k+l=m} choose(m,k) * (j/n^2) * B[i,k] * u1[j,l]`.

Output row `0` is exactly zero. Expanding the two mixed-scale factors and collecting only like theorem scales gives:

- ordinary `Lambda^0 ... Lambda^-4`;
- pressure-linear `a^2 Lambda^-1 ... Lambda^-3`;
- pressure-square `a^4 Lambda^-2`.

The `a^2` and `a^4` pieces stay in signed-log form. Only their finite Decimal normalized numerators are combined by the bilinear identity, so mathematically nonzero theorem-scale terms cannot silently underflow to binary64 zero.

## Validation

Regression coverage verifies:

1. The constructor accepts only the genuine theorem-scale first Picard state and keeps complete-`slow2`, `naturalRemainder(x1)`, and fixed-point flags false.
2. The ordinary `Lambda^0` component agrees exactly with a separate reference-only reconstruction of `average`, multiplication by `2*D*eta`, and the literal pinned `dot1` double sum.
3. Every ordinary, `a^2`, and `a^4` numerator agrees exactly with an independent literal mixed-scale `dot1` expansion.
4. `dot1` row zero is exactly zero while still exercising the genuine upstream state validation.
5. Actual nonzero pressure-linear terms retain common log scale `2*log(a)` with explicit inverse-`Lambda` powers, and the pressure-square term retains `4*log(a)` with explicit `Lambda^-2`.
6. Current nonzero theorem-scale pressure pieces fail closed on binary64 underflow rather than becoming zero.
7. Invalid indices, out-of-window eta, and surrogate/non-first-Picard inputs fail closed.

## Remaining boundary

This does not complete `slow2(x1)`. The remaining `average(u1)`-coupled branches through pinned `mixed1` and `param1` still require the same mixed-scale treatment, as do the remaining angular and source/pressure/resolvent operations needed for complete `naturalRemainder(x1)`. Only then can a genuine `x2` be materialized and the pinned contraction/convergence argument be attached to obtain fixed-point `phi/u`. Derived average/pressure, `NaturalProfileAssembly`, global all-index weighted `AxisSpace` certification, and support/moment/matching/cone closure remain downstream.
