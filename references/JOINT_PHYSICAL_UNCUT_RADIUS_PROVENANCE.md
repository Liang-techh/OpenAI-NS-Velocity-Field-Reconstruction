# Finite joint physical uncut-radius provenance (Agent 7 / Issue #2)

This increment is downstream of the provider-owned exact-majorant, common-support,
exact-recurrence, recursive SlowBorel/DiagonalScale schedule, joint ordinary-jet,
and joint physical-jet target-prefix chain. It does not implement or duplicate
Agent 2 hierarchy jets/coefficient recurrences.

## Source pin

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Modules:
  - `NavierStokes/GenericRealizationBounds.lean`
  - `NavierStokes/SmoothCutoffs.lean`
  - `NavierStokes/SlowBorelBase.lean`
- Pinned facts:
  - `GenericRealizationBounds.prefix_eventually_uncut`
  - `GenericRealizationBounds.realized_tail_bound`
  - `SmoothCutoffs.scaledCutoff_one_of_abs_le`
  - `SlowBorelBase.exists_physical_uncut_tail`

For a positive monotone cutoff schedule `a_j`, the pinned realization theorem
prints the exact finite-prefix plateau condition

`0 < q < 1/(2*a_J)`.

The already-landed recursive schedule satisfies the stronger doubling envelope
`a_{j+1} >= 2*a_j`, so every retained stage `j<=J` obeys `a_j<=a_J`. Therefore,
inside the strict open interval above, every retained cutoff argument satisfies
`|a_j q|<1/2` and the pinned smooth cutoff is exactly one.

## Increment

`background_joint_physical_uncut_radius.py` consumes the same provider-owned
physical target prefix selected by the parent stack. It does not accept a
caller radius or a second schedule. It extracts the actual recursively selected
integer scale `a_J` and records the theorem-ready radius exactly as

`delta_J = 1/(2*a_J)`.

The certificate checks positivity, monotonicity, the recursive doubling
envelope, exact endpoint arithmetic `a_J*delta_J=1/2`, and exact retained-stage
edges `a_j*delta_J<=1/2`. `certify_q` accepts only exact integer/Fraction values
and requires the strict theorem window before returning every exact `a_j*q`.
No cutoff function is sampled numerically.

Regression uses the parent physical target request `h=1/8`, `M=1`, `P=1/8`,
`Jmin=4`, `initial_lower_bound=3`. It selects `J=16`, obtains the actual
recursive scale `a_16=51010836299776`, and therefore the exact open plateau
radius `1/102021672599552`. At half that radius the final cutoff argument is
exactly `1/4`; equality at the radius itself is rejected because the pinned
tail theorem uses a strict neighborhood condition.

## Fail-closed truth boundary

This increment verifies the finite **plateau arithmetic and theorem hypotheses**
of the selected recursive prefix. Python does not replay Lean's eventual-equality
proof, materialize the infinite `slowSum`, prove the physical-chart finite
constant, evaluate the actual PDE residual, quantify over arbitrary hierarchy
orders, or prove Proposition 5.3.

Truth boundary:

- `finite_uncut_plateau_arithmetic_verified = true`
- `one_provider_owned_schedule_used = true`
- `function_level_uncut_germ_verified = false`
- `physical_chart_finite_bound_verified = false`
- `actual_slow_sum_materialized = false`
- `all_order_hierarchy_verified = false`
- `infinite_diagonal_schedule_verified = false`
- `pde_residual_tail_verified = false`
- `all_jets_flat = false`
- `super_algebraic = false`
- `paper_exact = false`
- `full_reconstruction = false`
- `paper_exact_velocity_available = false`

## Exact upstream return to Agent 2

The new exact plateau radius removes a downstream finite-cutoff ambiguity but
not the central hierarchy frontier. Agent 2 still must expose a real assembled
coefficient object, theorem-backed derivative majorants, common compact support,
and exact retained-recurrence evidence for every positive order requested by
arbitrary physical targets. The current finite PositiveAxis/Picard derivative
stack cannot be extrapolated into that total provider.
