# Section 8 exact five-row algebra provenance

Status: `formal-structure` only. This record does not promote any field to paper-exact status.

## Pinned formal source

- Repository: `https://github.com/openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/FiveRowRank.lean`
  - `NavierStokes.FiveRowRank.FiveRows`
  - `angularPowers`, `axialPowers`, `angularDebt`, `axialDebt`
- `NavierStokes/MeanRankUpdate.lean`
  - `NavierStokes.MeanRankUpdate.normalizeDebt`
  - `NavierStokes.MeanRankUpdate.physical_five_rows`

## Executable increment

`section8_exact_five_row_algebra.py` replays the exact finite-dimensional debt normalization and five-row scaling algebra with `fractions.Fraction`. It certifies, without tolerances, the two zero auxiliary moments and the three physical debt-cancelling targets `(0, 0, -P, -J_theta, -J_z)`. Float and Decimal inputs are rejected so numerical near-zero cannot substitute for an identity.

This complements, but does not replace, the existing diagnostic `mean_rank_update.py`: that module evaluates an explicit repository bump and numerically inverts its moment matrix, whereas the pinned Lean localized repair uses a noncomputable smooth bump. The new exact certificate therefore makes no claim that a localized radial inverse or compact correction field has been materialized.

## Truth boundary

The following remain false:

- `localized_moment_inverse_materialized`
- `actual_wave_defect_consumed`
- `actual_compact_correction_field_materialized`
- `paper_exact_velocity_available`
- Section 9 infinite correction-sequence convergence
- Eq. (9.21) infinite summed local field

The required upstream 3→8 handoff remains a source-revision-bound materialized Agent-3 oscillatory wave/defect/stage producer. Once that exists, its actual debt can be passed through the Section 8 localized repair rather than through a surrogate or caller-invented vector.
