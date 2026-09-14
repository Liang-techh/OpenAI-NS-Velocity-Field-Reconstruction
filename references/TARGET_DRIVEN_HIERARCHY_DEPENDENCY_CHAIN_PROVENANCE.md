# Target-driven hierarchy dependency-chain provenance

## Scope

This Agent-7 / Issue-#2 increment closes one provenance seam in the finite target-driven Section 5 closure path.  The existing target-driven prefix already obtains every finite `C[j,m]` row and every retained recurrence identity from one hierarchy provider, and the supported-prefix layer additionally requires the positive-order compact-support witness to share an own-order coefficient artifact with each derivative-bound row.  This increment requires one *single literal* own-order `(coefficient_order, artifact, provider)` key to occur simultaneously in all three evidence families at every selected positive order:

1. every `C[j,m]` row, `0 <= m <= j+2`;
2. the pinned common-support witness for the assembled positive-order coefficient; and
3. the exact retained recurrence identity at order `j`.

The purpose is to prevent a formally exact recurrence for a different coefficient object from being cross-wired into derivative/support evidence that merely carries the same hierarchy/source/coefficient-state labels.

## Pinned source

The formal source remains `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, especially the Section 5 assembled-coefficient support/compactness results, the SlowBorel template-jet bound and admissible recursive scales, the retained slow-expansion recurrence identities, and `DiagonalScale.exists_uniform_tail_order`.

The production path in this increment composes only already-existing repository certificates:

- `background_target_driven_hierarchy_prefix.py` for target-selected hierarchy-owned `C[j,m]` rows, the recursive SlowBorel/DiagonalScale schedule, and exact retained recurrence identities;
- `background_target_driven_support_prefix.py` for hierarchy-owned common support of the same finite positive-order prefix; and
- `background_recurrence_cancellation.py` for exact rational formal cancellation, where any nonzero rational defect remains nonzero.

No hierarchy jet, coefficient recurrence, or Agent-2 derivative constructor is reimplemented here.

## Fail-closed contract

`bind_target_driven_dependency_chain(...)` intersects the own-order dependency keys from support, the retained recurrence, and every derivative-bound row.  Construction fails if the intersection is empty.  Pairwise overlap is intentionally insufficient: for example, a support witness carrying artifacts `A` and `B`, `C[j,0]` carrying `A`, later derivative rows carrying `B`, and the recurrence carrying `A` does **not** establish one coefficient object behind all evidence and is rejected.

The one-shot `certify_target_driven_dependency_chain(...)` first invokes the existing provider-owned target-prefix constructor.  Callers therefore cannot supply a separate `C[j,m]` table, support table, or retained-recurrence list.  Provider order/support/recurrence frontiers propagate as failures and no partial certificate is returned.

Retained cancellation is still decided before this layer by exact `Fraction` algebra.  The regression injects a defect of exactly `10^-30`; it is rejected rather than treated as numerically zero.

## Truth boundary

This certificate applies to one finite prefix chosen for one exact target request.  It does not prove that the hierarchy provider is total for arbitrary coefficient order, does not construct an infinite diagonal schedule, and does not prove that the actual Navier--Stokes residual is controlled by the scalar tail.  In particular it does **not** establish Proposition 5.3, all-jets flatness, super-algebraic residual convergence, paper-exact velocity, or full reconstruction.

The corresponding truth flags remain:

- `finite_target_prefix_only = true`
- `one_owned_artifact_per_positive_order = true` within the selected finite prefix
- `retained_recurrences_exact = true` within that prefix
- `all_order_hierarchy_verified = false`
- `infinite_diagonal_schedule_verified = false`
- `pde_residual_tail_verified = false`
- `all_jets_flat = false`
- `super_algebraic = false`
- `paper_exact = false`
- `paper_exact_velocity_available = false`
- `full_reconstruction = false`

## Exact upstream return to Agent 2

The arbitrary-order hierarchy API needed by this downstream lane must eventually expose, for every requested positive order, the actual assembled coefficient object/artifact such that **all** hierarchy-owned analytic derivative bounds, its common compact-support provenance, and the exact retained recurrence decomposition name that same own-order artifact under one hierarchy/source/coefficient state.  A finite derivative jet frontier, or recurrence evidence tied only by labels rather than by the actual coefficient artifact, cannot discharge this contract.

Agent-2 PRs #319/#322 advance hierarchy-owned third-eta lower-source/PositiveAxis forcing composition, but they remain finite derivative-chain increments.  They do not yet expose a total arbitrary-order assembled-coefficient provider with the common dependency identity required here, so this increment deliberately remains finite/fail-closed rather than extrapolating those jets to all orders.
