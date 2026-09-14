# Finite joint-jet target-prefix provenance (Agent 7 / Issue #2)

This increment is downstream of the hierarchy-owned exact-majorant / support /
retained-recurrence chain. It does not implement Agent 2 hierarchy jets or
coefficient recurrences.

## Source pin

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Module: `NavierStokes/SlowBorelBase.lean`
- Relevant pinned facts:
  - `SlowBorelBase.exists_admissibleScales`
  - `SlowBorelBase.ordinary_tail_bound`
  - `SlowBorelBase.cutPrefix_eventually_uncut`
  - `SlowBorelBase.exists_ordinary_uncut_tail`
  - the underlying `DiagonalJetBounds.exists_prefix_gain`

The pinned ordinary-tail estimate is

`||D^m(slowSum-cutPrefix_J)|| <= 2^-J q^(h(J+1)-m)`

for `m <= J+3`. `exists_ordinary_uncut_tail` then chooses one `J` for an
entire finite derivative budget `m <= M` and target power `P`.

## Increment

`background_joint_jet_target_prefix.py` replays only the exact common-prefix
order arithmetic, on top of the existing provider-owned exact-majorant chain.

For exact positive `P`, finite `M`, and requested lower bound `Jmin`, it runs the
existing target-driven chain once at the weakest derivative exponent `m=M`,
while enforcing `J >= max(M,Jmin)`. Therefore every requested `m <= M` has

`h(J+1)-m >= h(J+1)-M >= P`.

The selected prefix is one shared recursive SlowBorel/DiagonalScale schedule;
it is not a collection of independently selected per-derivative schedules. All
`C[j,m]` rows remain hierarchy-provider-owned exact rational majorants rounded
outward before the executable scheduler, and the parent chain still requires
the same actual coefficient artifact to own derivative bounds, common compact
support, and exact retained recurrence evidence.

## Fail-closed boundary

The constructor does **not** replay the full function-level
`exists_ordinary_uncut_tail` theorem. In particular it does not:

- materialize an infinite hierarchy or `slowSum`;
- prove the local `cutPrefix_eventually_uncut` germ or construct its `delta`;
- prove a runtime PDE residual estimate;
- quantify over every `M` and `P`;
- infer all-jets-flatness or super-algebraic convergence from one finite request.

A provider frontier anywhere before the selected `J` aborts the whole joint
request; easier lower derivative orders are not returned as a partial success.
Float target powers are rejected before the theorem gate.

## Regression

`tests/test_background_joint_jet_target_prefix.py` verifies a nontrivial common
prefix with `h=1/8`, `M=1`, `P=1/8`, where the exact minimal target index above
`Jmin=4` is `J=8`. The two exact ordinary powers are `9/8` and `1/8`, with
the common dyadic factor `2^-8`. Separate regressions check `Jmin` domination,
provider-frontier failure, float-target rejection, and invalid finite budgets.

## Truth boundary

- `joint_finite_jet_target_arithmetic_verified = true`
- `one_common_recursive_schedule_for_requested_jets = true`
- `finite_joint_jet_prefix_only = true`
- `uncut_tail_germ_verified = false`
- `all_order_hierarchy_verified = false`
- `infinite_diagonal_schedule_verified = false`
- `pde_residual_tail_verified = false`
- `all_jets_flat = false`
- `super_algebraic = false`
- `paper_exact = false`
- `full_reconstruction = false`
- `paper_exact_velocity_available = false`

## Exact upstream return to Agent 2

The next real upstream requirement is unchanged but now has the finite
quantifier shape made explicit: for every requested positive coefficient order
needed by the common `J(M,P)`, Agent 2's eventual assembled-coefficient API must
materialize the actual coefficient object and theorem-backed derivative
majorants/support/exact recurrence evidence under one coefficient state. The
current third-eta / first-Picard frontier (#319/#322/#326) cannot yet answer
arbitrary `J(M,P)`, so this certificate intentionally fails closed beyond that
frontier.
