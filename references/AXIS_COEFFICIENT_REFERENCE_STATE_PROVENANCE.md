# Actual reference coefficient-state provenance

Status: **formal-structure**. This increment does **not** make the leading profile or velocity paper-exact.

## Pinned official source

Official Lean repository: `openai/NavierStokesAndEuler` at commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The executable adapter in `src/openai_ns_reconstruction/axis_coefficient_reference_state.py` follows these pinned definitions:

- `NavierStokes/AxisCoefficientSpace.lean`
  - `RawJets`: coordinates indexed by `(radial degree n, parameter derivative m, eta)`.
  - `jet`: the stored normalized coordinate multiplied by `AxisWeightEstimates.weight`.
  - `Compatible`: adjacent actual parameter jets obey the FTC identity throughout the compact window.
  - `ofJetFamily`: a continuous family with actual adjacent derivatives and a global weighted bound gives a genuine element of the complete coefficient space.
- `NavierStokes/NaturalAxisCoefficients.lean`
  - `window = [-11/10, 11/10]`.
- `NavierStokes/AxisContraction.lean`
  - `referencePair` is the pair about which the nonlinear contraction is run.

## Landed executable inputs

This adapter does not invent coefficient values. It combines the already-landed actual-schedule derivative families:

- `axis_reference_angular_jets.py`: arbitrary finite actual eta-jets of the angular reference coefficients, derived from the certified schedule sigma and the pinned radial recurrence.
- `axis_reference_axial_jets.py`: arbitrary finite actual eta-jets of `u0[1]=-(1/2) inverseL*zStar`, including the actual SchedulePressure derivative chain.

Both paths derive `epsilon=rho/2` from the same actual-schedule analytic-neighborhood certificate. The adapter checks that the angular/axial epsilon, sigma, schedule data, and `j` agree before exposing the paired state.

`AxisCoefficientJetState.jet(n,m,eta)` is the **unnormalized actual derivative**, while `normalized_coordinate` divides by the exact Python realization of the pinned axis weight. No caller-supplied coefficient table, sigma, pressure, rho, or epsilon is accepted by `actual_schedule_reference_axis_state`.

## Independent regression boundary

`tests/test_axis_coefficient_reference_state.py` checks the two component families through an independent SciPy adaptive-quadrature FTC identity, verifies the zeroth jets against the previously landed reference-pair coefficient evaluator, verifies the pinned weight normalization directly from its formula, and checks the exact enlarged window guards.

The FTC regression is finite numerical evidence for the executable derivative chain; it is not a replacement for the Lean proof.

## What remains unresolved

This increment deliberately stops before claiming a complete `AxisCoefficientSpace.AxisSpace` backend. In particular it does **not** yet provide:

1. a global all-index weighted supremum norm certificate needed to instantiate the full Banach-space membership computationally;
2. executable `AxisOperators` / `coefficientOperators` on the state representation;
3. the official `naturalRemainder` evaluation at the actual reference pair;
4. the first genuine Picard iterate or the fixed-point fields `phi/u`;
5. the derived average/pressure fields or `NaturalProfileAssembly`;
6. support, moments, matching, and cone-condition closure for Theorem 4.6.

Accordingly `paper_exact` and `paper_exact_velocity_available` remain fail-closed (`false`).
