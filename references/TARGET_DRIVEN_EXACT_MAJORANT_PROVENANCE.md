# Target-driven exact-majorant provenance (Agent 7 / Issue #2)

This increment sits downstream of the existing Agent-7 finite target-driven chain and
upstream of any all-order claim.  It does **not** implement Agent 2 hierarchy jets
or coefficient recurrences.

## Source pin

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Relevant modules:
  - `NavierStokes/SlowBorelBase.lean`
  - `NavierStokes/DiagonalScale.lean`
  - `NavierStokes/AssembledSlowBase.lean`
  - `NavierStokes/SlowExpansionResidual.lean`
- Relevant pinned facts already replayed by the parent stack:
  - `SlowBorelBase.exists_template_jet_bound`
  - `SlowBorelBase.exists_admissibleScales`
  - `DiagonalScale.exists_uniform_tail_order`
  - `AssembledSlowBase.extendedCoefficient_support`
  - `AssembledSlowBase.extendedCoefficient_compactSupport`

The formal source supplies analytic upper-bound existence and the recursive cutoff
construction.  It does not say that its existential real constants are rational.
The exact `Fraction` in this increment is therefore an implementation-level
**certified rational upper majorant** that the eventual real hierarchy backend must
justify for the same actual coefficient artifact.

## Closed seam

The executable SlowBorel scheduler currently carries `C[j,m]` as binary64 floats.
A nearest-float conversion can round a valid positive upper majorant downward
*before* the recursive scale inequality is evaluated.  This increment forbids that
path:

1. the hierarchy provider owns `ExactHierarchyJetMajorant(j,m)` as an exact positive
   integer/Fraction, plus hierarchy/source/coefficient-state and dependency provenance;
2. binary64 materialization is outward: if nearest conversion lies below the exact
   rational, it advances with `nextafter(..., +inf)`;
3. the converted value is rechecked exactly with `Fraction.from_float`;
4. non-finite/overflowing conversion fails closed;
5. the existing target-driven support/recurrence dependency chain is then executed
   unchanged, so the same own-order coefficient artifact must still own every bound,
   support witness, and exact retained recurrence row.

No caller-supplied `C[j,m]` table is introduced.

## Regression

`tests/test_background_target_driven_exact_majorants.py` includes a rational
`1 + 2^-60`, whose nearest binary64 value is exactly `1.0`; the production
conversion must step to `nextafter(1.0, +inf)` and prove that the resulting float
dominates the rational.  The suite also rejects floating "exact" inputs, mismatched
coefficient state, binary64 overflow, and a target that requests an order beyond a
real provider frontier.

## Truth boundary

This is one finite target-selected arithmetic/provenance certificate only:

- `finite_target_prefix_only = true`
- `outward_binary64_majorants_verified = true`
- `all_order_hierarchy_verified = false`
- `infinite_diagonal_schedule_verified = false`
- `pde_residual_tail_verified = false`
- `all_jets_flat = false`
- `super_algebraic = false`
- `paper_exact = false`
- `full_reconstruction = false`
- `paper_exact_velocity_available = false`

## Exact upstream return to Agent 2

For every requested positive coefficient order, the eventual hierarchy-owned API
must produce a theorem-backed upper majorant for each required `C[j,m]` that is
bound to the same actual assembled coefficient artifact as its support and exact
retained recurrence evidence.  A raw binary64 estimate is insufficient for this
path because it cannot certify that rounding preserved the upper-bound direction.
Finite eta-derivative frontier PRs remain useful inputs, but they are not a total
arbitrary-order provider and therefore cannot trigger an all-order certificate.
