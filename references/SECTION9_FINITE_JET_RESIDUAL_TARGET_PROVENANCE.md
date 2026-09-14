# Section 9 finite-jet residual target provenance

Status: **formal-structure / exact finite arithmetic only**.

## Pinned source

- Formal repository: `https://github.com/openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Primary file: `NavierStokes/ActualIterationLedger.lean`
- Symbols:
  - `NavierStokes.ActualIterationLedger.residualRate`
  - `NavierStokes.ActualIterationLedger.residualRate_tendsto_atTop`
  - `NavierStokes.ActualIterationLedger.eventually_residualRate_ge`
- Supporting exact formulas are replayed by the stacked `section9_residual_rate_target.py` path from `CoordinateAlgebra.A`, `PhysicalGraphBounds.graphLoss`, `ActualIterationLedger.residualLoss`, and `ActualIterationLedger.residualWave`.

## Executable increment

`section9_finite_jet_residual_target.py` takes exact `h`, `beta`, a finite derivative budget `M`, and one exact target residual-rate exponent.  It obtains the individually minimal target-reaching stage for every `m=0,...,M`, chooses one common stage equal to their maximum, and then rechecks all derivative rows at that common stage.

Global minimality is checked independently: when the selected common stage is positive, at least one requested derivative row must still miss the target at the immediately preceding stage.  The implementation does **not** assume derivative-order monotonicity; it inspects every finite row, so no hidden sign assumption on `beta` is introduced.  Exact `Fraction` arithmetic is retained throughout and float/Decimal theorem inputs are rejected.

This is a finite `C^M` exponent-budget constructor.  It is useful downstream because an actual Section 9 stage producer can eventually be asked for one prefix deep enough to meet a whole finite jet budget instead of choosing a different prefix for every derivative order.

## Truth boundary

The repository still does not have an Agent-3 source-revision-bound materialized oscillatory wave/defect/stage producer.  Therefore this increment does **not** consume an actual stage field and does not evaluate an actual Navier--Stokes residual.

The following remain false:

- `actual_stage_field_consumed`
- `actual_residual_evaluated`
- `all_derivative_orders_certified`
- `infinite_correction_sequence_certified`
- `eq_9_21_summed_field_certified`
- `paper_exact_velocity_available`
- `full_reconstruction`

A finite derivative budget, even for an arbitrarily large chosen `M`, is not an all-jets statement and is not a convergence proof for Eq. (9.21).  The next upstream requirement from Agent 3 remains an actual wave/defect/stage producer tied to the same paper source revision; only then can these exact target indices be checked against genuine stage residuals and fed toward an Agent-4 local field/residual contract.
