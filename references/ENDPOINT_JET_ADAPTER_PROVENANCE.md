# Section 10 full-spacetime to normal-jet adapter provenance

Status: **formal-structure**. This increment connects the full residual-derivative limit interface used by the pinned `CandidateFromLimits` theorem to the normal time-jet interface consumed by the repository's Taylor--Borel evaluator. It does not construct the missing limits or certify endpoint smoothness.

## Official source mapping

Pinned official Lean repository:

`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant modules and statements:

- `NavierStokes/CandidateFromLimits.lean`: the hypothesis `hlim` supplies, for every order `n`, locally uniform left limits `L x n` of the **full spacetime iterated Frechet derivative** of the actual Navier--Stokes residual.
- `NavierStokes/SpacetimeGluing.lean`: `timeVector : SpaceTime := (1, 0)` and `normalIter` repeatedly differentiate in that time direction. `smoothExtension` feeds the resulting `normalTrace` into `SpatialBorelExtension.rightExtension`.
- `NavierStokes/SpatialBorelExtension.lean`: the right extension accepts one spatially smooth normal time-jet function for each Taylor degree and chooses degree-dependent scales from derivative/localization bounds.

Thus the Borel coefficients used at degree `n` are not arbitrary components of `L x n`: they are the repeated contraction of the full order-`n` tensor with the time direction in every derivative slot.

## Executable representation

`src/openai_ns_reconstruction/endpoint_jets.py` uses the dense convention

`shape = (3,) + (4,)*n`,

where the first axis is the three-vector output and every derivative slot is ordered `(time, x, y, z)`. Under this convention the repeated contraction with Lean's `timeVector=(1,0)` is exactly

`J[:, 0, ..., 0]`.

`normal_time_jet_from_full` performs and validates this contraction. `normal_jet_family` adapts a full-jet family to the existing `BorelRightExtension` API, and `borel_extension_from_full_spacetime_jets` composes the two interfaces without inventing any scale or jet data.

## Independent check

`tests/test_endpoint_jets.py` uses the analytic spacetime polynomial

`F(t,x,y,z) = (t+x, t^2+y, t^3+z)`

at `t=1`. It populates full derivative tensors, including nonzero spatial entries that must not leak into the normal contraction, verifies the expected first three normal time derivatives, then passes those full tensors through the adapter and existing Borel evaluator. At `t=1.04`, all nonzero polynomial Taylor terms lie on the exact cutoff plateau, so the resulting future branch is checked directly against the analytic polynomial rather than against a second implementation of the production sum. Shape, non-finite and invalid-degree inputs fail closed.

## Truth boundary

This adapter does **not** establish any of the analytic hypotheses required by `CandidateFromLimits` or `SpacetimeGluing`. In particular it does not:

- show that supplied tensors are derivatives of the actual localized Navier--Stokes residual;
- prove `TendstoLocallyUniformly` as `t -> 1-` for any derivative order;
- derive the full derivative recurrence or closed-past joint smoothness;
- construct the `SpatialBorelExtension.localScale` derivative-majorant witnesses;
- prove that the repository's explicit cutoff transition agrees pointwise with Mathlib's noncomputable `ContDiffBump`;
- provide the still-missing completed local velocity/pressure from Issues #1--#3.

Therefore Stage 7 remains **formal-structure**, and `paper_exact_velocity_available` must remain `false`. The next substantive endpoint step is to obtain actual residual derivative limits (starting with a reusable numerical/analytic certificate interface if the upstream field is still unavailable), not to populate this adapter with arbitrary tensors and call them paper-exact.
