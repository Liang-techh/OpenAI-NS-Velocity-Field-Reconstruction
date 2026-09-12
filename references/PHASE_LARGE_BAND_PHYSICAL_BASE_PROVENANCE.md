# Section 6→7 large-band physical-base identity provenance

Status: **formal-structure only**. This artifact does not materialize a paper
background field and does not make the paper-exact velocity available.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- File: `NavierStokes/PrimaryGeometryAssembly.lean`
- Field-level declarations pinned by the adapter:
  - `PrimaryGeometryAssembly.construction_frequency`
  - `PrimaryGeometryAssembly.construction_axial`
  - `PrimaryGeometryAssembly.carrier_time_positive`
  - `PrimaryGeometryAssembly.Prepared.radius_pos`
  - `PrimaryGeometryAssembly.frequency_eq_physical`
  - `PrimaryGeometryAssembly.axial_eq_physical`
  - `FinalSlowBase.velocity`
  - `BaseChartJets.cellBand`

For a `Prepared` object `a` and active index `L`, the formal source identifies
`construction ... c` frequency/axial fields with `frequency`/`axial` at
`cellBand L`, independently for either `c : Fin 2`. On the active carrier,
`carrier_time_positive` and `a.radius_pos` supply the positivity hypotheses used
by the physical identities. `frequency_eq_physical` and `axial_eq_physical`
then identify those normalized fields with scaled components 1 and 2 of the
same `FinalSlowBase.velocity` evaluated at the corresponding band point.

## Implemented claim

`phase_large_band_physical_base.py` adds a fail-closed theorem-witness bridge
from an existing `LargeBandPreparedSourceBinding` to that exact field identity
chain. Admission requires exact box, source revision, and opaque Prepared
instance identity, plus separate theorem certifications for:

1. `Prepared.base` using the same frequency/axial field family;
2. `construction_frequency` and `construction_axial`;
3. carrier time and radius positivity needed on the active domain;
4. `frequency_eq_physical` and `axial_eq_physical`; and
5. reuse of the same physical base by both sign branches.

Only `analytic-theorem` and `formal-theorem` evidence is accepted. Sampled,
fitted, and numeric-scan evidence is rejected. Tests use synthetic theorem
metadata fixtures only; they are not manuscript field values or an exported
Lean application.

## Explicit non-claims

The Python adapter does not evaluate `FinalSlowBase.velocity`, `bandPoint`,
`Prepared.base`, `frequency`, or `axial`; it does not independently discharge
the theorem hypotheses; and it does not infer field equality from sampled
values. Therefore a successful admission means only that an external theorem
application has been fail-closed against the pinned identity chain.

The following remain false:

- `actual_physical_base_values_materialized`
- `actual_base_fields_verified`
- `uniform_eq_7_9_to_7_11_verified`
- `paper_exact_velocity_available`

## Remaining paper-exact boundary

A machine-linked export must still provide the concrete `Prepared` theorem
application and actual active `CellIndex` correspondence, so these certified
identities can be instantiated for every active carrier rather than represented
by theorem metadata. Those field-level witnesses must then feed the existing
phase/base/frame/damping coverage gates for Eqs. (7.9)–(7.11). Stress-cone
amplitudes, curl-realized divergence-free waves, compact mean corrections, and
Section 9 residual improvement remain separate later gates.
