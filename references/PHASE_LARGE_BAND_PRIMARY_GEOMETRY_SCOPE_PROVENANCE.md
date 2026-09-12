# Section 7 large-band primary-geometry scope provenance

Status: **formal-structure**.

This increment strengthens the upstream side of the finite active-family
coverage gate.  `phase_large_band_family_coverage.py` deliberately accepts an
externally declared exact finite signed family scope; it checks complete
coverage of that declaration, but it does not identify the declaration with
the formal primary-geometry index used by the reconstruction.

The pinned formal source at `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd` provides the relevant structural
objects:

- `BaseChartJets.CellIndex` in `NavierStokes/BaseChartJets.lean` is defined as
  the subtype of actual positive active labels whose band is at least the
  fixed lower band `N`.
- `PrimaryGeometryAssembly.cellDomain` in
  `NavierStokes/PrimaryGeometryAssembly.lean` uses that index, with scale
  `ChartScales.S (BaseChartJets.cellBand L)` and carrier given by the formal
  open cell.
- `PrimaryGeometryAssembly.family` constructs
  `BasePhaseGeometry.FamilyData (domain W a.N) ...` for `c : Fin 2`, so the
  two phase-family branches share the same underlying primary-geometry index.

`phase_large_band_primary_geometry_scope.py` does **not** attempt to evaluate
or enumerate those noncomputable Lean objects in Python.  Instead it accepts a
theorem-provenanced export of one fixed-band unsigned box slice.  The witness
must pin all three definitions above, certify that the supplied box manifest
is complete for that band, certify the band/cell projection and carrier-domain
identity, certify that `PrimaryGeometryAssembly.family` uses the same index,
and certify a bijective identification of its two `Fin 2` branches with the two
Section 6 signs.  No particular `Fin 2 -> sigma` ordering is hard-coded.

The Python verifier independently performs the deterministic set check:
for every unsigned box `a` it reconstructs exactly
`(ell,a,-1)` and `(ell,a,+1)` and requires the resulting set to equal the
existing `LargeBandActiveFamilyScopeWitness` exactly.  Missing boxes, extra
boxes, duplicate boxes, band/source drift, stale Lean pins, or sampled/fitted/
numeric-scan evidence all fail closed.

The tests use a synthetic theorem metadata fixture and explicitly do **not**
claim that its two boxes are manuscript data or a real `CellIndex` export.
Therefore successful admission only establishes consistency between an
external theorem witness and the declared Python scope.  It does not itself
prove or extract the paper's actual active set.

Truth boundary remains:

- `actual_cell_index_enumeration_machine_verified = false`
- `actual_active_family_scope_machine_verified = false`
- `actual_base_fields_verified = false`
- `uniform_eq_7_9_to_7_11_verified = false`
- `paper_exact_velocity_available = false`

No surrogate background, sampled uniformity claim, amplitude, or oscillatory
wave is introduced.
