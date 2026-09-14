# Section 8 compact support geometry provenance

## Scope

This increment machine-replays only the support geometry used by the Section 8 five-row compact mean correction. It does **not** materialize the smooth correction functions, solve their generalized-power moment matrix numerically, consume an oscillatory-wave defect, or promote the reconstruction to paper-exact status.

The production certificate works in the exact affine coordinate

`R = (1 - alpha) * a + alpha * b`

instead of evaluating the possibly noncomputable real patch endpoints `a,b`. Consequently all cell/support ordering checks use `fractions.Fraction`; binary floats and caller-selected generic cutoff geometry are rejected.

## Pinned formal source

Repository: `https://github.com/openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Files:

- `NavierStokes/FiveRowRank.lean`
- `NavierStokes/LocalizedMomentRepair.lean`

Replayed definitions/theorems:

- `NavierStokes.FiveRowRank.cellStep`
- `NavierStokes.FiveRowRank.cellLower`
- `NavierStokes.FiveRowRank.cellUpper`
- `NavierStokes.FiveRowRank.cell_separated`
- `NavierStokes.FiveRowRank.cell_union_subset`
- `NavierStokes.LocalizedMomentRepair.innerLower`
- `NavierStokes.LocalizedMomentRepair.innerUpper`
- `NavierStokes.LocalizedMomentRepair.repair_tsupport_subset_open`

For a family of `n` correction cells the exact affine coordinates are

- lower: `(2*j + 1)/(2*n + 1)`
- upper: `(2*j + 2)/(2*n + 1)`

and the actual localized repair is supported in the middle half of each cell. Section 8 fixes `n=3` for the angular repair and `n=2` for the axial repair.

The certificate checks, with exact rational arithmetic, that every cell lies strictly inside `(a,b)`, every middle-half support lies strictly inside its cell, and adjacent cells retain the positive unused gap required by the generalized-power moment inverse. It rejects altered cells, altered middle-half supports, floating-point coordinates, wrong 3/2 family sizes, and altered pinned provenance.

## Why this is not another numerical five-row solve

`src/openai_ns_reconstruction/mean_rank_update.py` remains a diagnostic floating-quadrature representative. The formal source instead constructs `LocalizedMomentRepair.repair` from an exactly nonsingular generalized-power moment matrix and proves its moment identities in Lean. This increment deliberately does not approximate those moment integrals and does not turn numerical near-zero values into identities.

The value of this increment is narrower: support preservation is now definition-exact and independent of the numerical diagnostic bump. A future materialized rank increment can be required to carry this fixed geometry rather than an arbitrary caller-supplied cutoff.

## Truth boundary and next blocker

Certified here:

- pinned three-angular/two-axial cell geometry;
- pinned middle-half compact support geometry;
- exact pairwise separation and strict containment as affine identities;
- exact dependency chain to the two pinned formal files.

Still **not** materialized here:

- actual reserved-patch endpoint values from the current wave/rank application;
- the Lean smooth cutoff values or generalized-power moment matrix inverse;
- exact five moment values as executable data;
- actual `rankIncrementState` / `rankStageState` / nonlinear remainder values;
- Agent 3's actual oscillatory-wave defect/stage fields;
- Section 9 infinite correction sequence or Eq. (9.21) limit field.

The precise 3 -> 8 handoff remains: Agent 3 must export the machine-materialized oscillatory wave/defect (or canonical stage producer) under the same source revision. Section 8 then needs a value-level producer for the actual rank increment/remainders, with the formal `repair_exact`, zero-mass, five-row, and support theorems attached. Only then can Section 9 consume actual stage fields rather than theorem identities.
