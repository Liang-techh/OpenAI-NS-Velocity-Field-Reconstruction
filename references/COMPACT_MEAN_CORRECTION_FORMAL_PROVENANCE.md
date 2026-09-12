# Compact mean correction formal provenance

Status: **formal-structure only**.

This layer is the paper-exact formal counterpart to the executable diagnostics in
`mean_rank_update.py`.  It does not reuse that module's compact polynomial bump,
which is intentionally a numerical surrogate rather than the noncomputable bump
in the pinned Lean development.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- File: `NavierStokes/CorrectionState.lean`
- Main theorem: `NavierStokes.CorrectionState.rank_rows_on_patch`
- Five-row dependency: `NavierStokes.MeanRankUpdate.physical_five_rows`

The formal source describes `rank_rows_on_patch` as the theorem that transfers
the compact five-row construction to the **actual base** when that base agrees
with the prescribed angular/axial background on the support of the constructed
bumps.  Its inputs include the actual `CorrectionState.debt c u n x` through the
rank desired profiles, together with positivity/nondegeneracy hypotheses for the
rank geometry.

## What the Python gate checks

`compact_mean_correction_formal.py` accepts only a `lean-formal-export` witness
for one exact theorem application.  It pins the repository, commit, file,
theorem and five-row dependency; records stable opaque identities for the rank
data, context, correction state, exact state debt, slow/plane points and the two
background models; and requires formal certification of the theorem's geometric
nondegeneracy and support-agreement hypotheses.

The gate rejects sampled, fitted and numeric-scan evidence.  It never evaluates
or reconstructs Lean fields in Python.

## Deliberate non-claim

The finite-head identity admitted by `actual_signed_mean_defect.py` is **not**
automatically identified with `CorrectionState.debt`.  No pinned theorem tying
those two objects has been admitted here, so this layer does not feed a synthetic
or inferred finite-head debt into the compact rank construction.

Accordingly all of the following remain false:

- `theorem_application_machine_replayed`
- `actual_debt_values_materialized`
- `actual_bump_values_materialized`
- `finite_head_signed_defect_link_verified`
- `compact_mean_correction_field_materialized`
- `paper_exact_velocity_available`

The regression fixture contains synthetic theorem metadata only.  It is not an
evaluated paper field, a reconstructed bump, or a replayed Lean application.
