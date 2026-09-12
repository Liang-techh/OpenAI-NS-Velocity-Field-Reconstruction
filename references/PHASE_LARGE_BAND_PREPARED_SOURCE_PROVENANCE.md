# Section 6→7 large-band `Prepared` source provenance

Status: **formal-structure only**.  This artifact does not make the paper-exact
velocity available and does not materialize the noncomputable Lean fields.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- File: `NavierStokes/PrimaryGeometryAssembly.lean`
- Declarations pinned by the adapter:
  - `PrimaryGeometryAssembly.Prepared`
  - `PrimaryGeometryAssembly.exists_prepared`
  - `PrimaryGeometryAssembly.prepared`
  - `PrimaryGeometryAssembly.family`
  - `PrimaryGeometryAssembly.construction`
  - `PrimaryGeometryAssembly.phases`

The formal source is materially stronger than a free source label.  `Prepared`
contains the LocalBase family plus polynomial frequency/axial jets and the
shared global constants.  `exists_prepared` constructs those data from the
actual `FinalSlowBase`/`BaseChartJets` chain and true-cone hypotheses;
`prepared` is `Classical.choice` of that theorem.  `family` sets its `base`
field to `a.base` and introduces the sign only through `phaseSign c`.
`phases` then applies `construction` to the same selected `prepared` object for
both `c : Fin 2` branches.

## Implemented claim

`phase_large_band_prepared_source.py` adds a fail-closed theorem-witness bridge
from an existing `LargeBandBaseSourceBinding` to that exact pinned constructor
chain.  Admission requires:

1. exact `(source_id, source_revision)` agreement with the already admitted
   sign-free slow-box base source;
2. a stable opaque `prepared_instance_id` supplied by the theorem/export layer;
3. theorem certification that the concrete `exists_prepared` application is
   valid and that the selected `prepared` object has the claimed identity;
4. theorem certification that `family ... .base` is the same source and that
   `phases` reuses the same prepared construction for both signs; and
5. exact repository/commit/file/declaration pins above.

Only `analytic-theorem` and `formal-theorem` evidence is accepted.  `sampled`,
`fitted`, and `numeric-scan` evidence is rejected rather than being upgraded to
a structural theorem witness.

The regression uses **synthetic theorem metadata fixtures only**.  Those
fixtures test the gate and are not manuscript data, a `Prepared` export, or a
paper background field.

## Explicit non-claims

The Python module does **not** discharge the hypotheses of `exists_prepared`,
evaluate `Classical.choice`, enumerate `BaseChartJets.CellIndex`, or materialize
`Prepared.base`, its jets, phase vectors, or oscillatory waves.  A caller's
certification flags are external theorem provenance; they are not converted
into an independent Python proof of the Lean application.

Accordingly the adapter keeps all of the following false:

- `actual_prepared_application_machine_verified`
- `actual_base_fields_verified`
- `uniform_eq_7_9_to_7_11_verified`
- `paper_exact_velocity_available`

## Remaining paper-exact boundary

A machine-linked export/theorem application must still bind concrete manuscript
parameters and hypotheses to `PrimaryGeometryAssembly.prepared` (or an
extensionally identified equivalent), connect the actual active
`BaseChartJets.CellIndex` family and concrete base/jet values, and feed those
field-level witnesses through the existing phase/base/frame/damping admissions.
Only after that can Eqs. (7.9)–(7.11) be upgraded from formal-structure
provenance toward an actual all-active-family verification.  Stress-cone
amplitudes, curl-realized divergence-free waves, compact mean corrections, the
Section 9 residual iteration, and the final paper-exact velocity remain separate
later gates.
