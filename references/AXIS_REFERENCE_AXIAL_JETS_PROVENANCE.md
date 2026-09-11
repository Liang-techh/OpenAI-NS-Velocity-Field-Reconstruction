# Actual Schedule Axial Reference Jets provenance

Status: **formal-structure**. `paper_exact_velocity_available=false`.

This increment advances the Stage-1 `AxisContraction.referencePair` from a zeroth-parameter-jet axial coefficient to arbitrary finite real-parameter derivative jets for the **actual constructed SchedulePressure datum**. It does not materialize the nonlinear fixed point or the final Theorem 4.6 profiles.

## Pinned source mapping

Official Lean repository/commit: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

- `NavierStokes/AxisContraction.lean`: `referencePair`, whose axial component is `-(1/2) J_1(inverseL * zStar)`.
- `NavierStokes/AxisWeightEstimates.lean`: `regularInverseJet`/`radialDivisor`, implying a radially constant datum contributes only radial degree one for `J_1`.
- `NavierStokes/NaturalAxisData.lean`: definitions of `L`, `U`, `H`, `d`, and `Z`.
- `NavierStokes/SchedulePressure.lean`: actual `axisPressure`, `axisPressure_contDiff`, `axisPressure_hasDerivAt`, and the actual-schedule natural-profile existence chain.
- `NavierStokes/PressureDatum.lean`: `kernel(a,eta)=(1+eta^2)^(-2a)` and the holomorphic/`ContDiff` pressure construction.
- `NavierStokes/AxisCoefficientSpace.lean`: stored parameter jets are genuine successive derivatives, not independent caller arrays.

## Executable realization

`src/openai_ns_reconstruction/schedule_axis_pressure_jets.py` computes normalized finite Taylor jets of the exact real pressure kernel by truncated analytic power-series algebra (`log` followed by `exp`). The pressure integral uses the same actual `OutgoingTail.finalAngular` clock and the official shape-exponent geometry: exponent one before flattening, variable exponent only on the compact flattening interval, and exponent zero after flattening. Constant-exponent regions reuse the landed plateau/tail reductions; only the compact flattening interval uses cached Gauss--Legendre quadrature.

`src/openai_ns_reconstruction/axis_reference_axial_jets.py` then differentiates the actual

`Z = -A(1-2 eta U)U - 4H - d P' + 4A eta P`

and `inverseL=1/L` by exact truncated Taylor algebra and constructs

`u0[1] = -(1/2) inverseL * zStar`,

with every radial degree other than one identically zero. The coefficient-space normalization uses the already landed official axis weight and the same actual-schedule analytic-neighborhood `epsilon=rho/2` certificate as the angular reference jets.

## Independent checks

Regression coverage checks the kernel Taylor series against closed-form special cases, checks pressure value/first derivative against the pre-existing independent `axis_pressure` / `axis_pressure_derivative` path, checks the axial zeroth jet against the landed `ActualScheduleReferencePair`, and checks first/second axial derivatives against symmetric finite-difference oracles of that separately implemented coefficient evaluator.

## Numerical and truth boundary

The new derivative chain is analytic in eta; it is **not** obtained from sampled polynomial fitting or finite-difference production code. However the compact flattening integrals are still floating-point Gauss--Legendre evaluations, just as the landed actual schedule pressure evaluator is. No rigorous quadrature enclosure or local Lean build certificate is claimed.

This increment therefore does **not** establish a complete compatible `AxisCoefficientSpace` state, the global weighted norm witness, `coefficientOperators`, `naturalRemainder`, a genuine Picard iterate/fixed point, `phi/u/average/pressure`, `NaturalProfileAssembly`, or final support/moment/matching/cone conditions. Stage 1 remains **formal-structure**, and `paper_exact_velocity_available` remains **false**.
