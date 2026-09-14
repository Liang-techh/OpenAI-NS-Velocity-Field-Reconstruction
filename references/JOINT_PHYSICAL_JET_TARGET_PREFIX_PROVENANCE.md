# Finite joint physical-jet target-prefix provenance (Agent 7 / Issue #2)

This increment is downstream of the hierarchy-owned exact-majorant, common-support,
exact-recurrence, and joint ordinary-jet prefix chain. It does not implement or
duplicate Agent 2 hierarchy jets/coefficient recurrences.

## Source pin

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Module: `NavierStokes/SlowBorelBase.lean`
- Pinned facts:
  - `SlowBorelBase.exists_ordinary_uncut_tail`
  - `SlowBorelBase.exists_physical_uncut_tail`
  - `SlowBorelBase.physicalChart_finite_bound`

The pinned physical theorem chooses the ordinary tail with requested power
`P + m`, then absorbs the physical-chart derivative loss. For one common finite
physical derivative budget `0 <= m <= M`, the worst exact scalar request is
therefore the ordinary request at `m=M` with target `P+M`.

The resulting exact exponent bookkeeping is

`ordinary exponent(m) = h(J+1)-m`

and, after spending the physical-chart derivative loss,

`physical exponent(m) = h(J+1)-2m >= P`.

## Increment

`background_joint_physical_jet_target_prefix.py` runs the already provider-owned
joint ordinary prefix once at target `P+M`. The same recursively selected
SlowBorel/DiagonalScale prefix, exact hierarchy majorants, common support, and
retained recurrence identities are shared by the whole finite physical jet
budget. The theorem-side requirement `0<h<1/2` is checked exactly through the
underlying rationalized `h`.

Regression uses `h=1/8`, `M=1`, `P=1/8`. The common minimal prefix is `J=16`;
the exact ordinary exponents are `(17/8, 9/8)` and the physical exponents after
the additional derivative loss are `(17/8, 1/8)`, with common dyadic factor
`2^-16`.

## Fail-closed boundary

This certificate does **not** construct or certify the physical-chart finite
constant `D`, the local uncut germ/delta, the infinite `slowSum`, the actual PDE
residual tail, an infinite diagonal schedule, or all-jets-flat/super-algebraic
convergence. A hierarchy frontier before the selected `J` aborts the whole
request. Binary-float target powers are rejected.

Truth boundary:

- `joint_finite_physical_jet_target_arithmetic_verified = true`
- `one_common_recursive_schedule_for_requested_physical_jets = true`
- `retained_recurrences_exact = true`
- `finite_physical_jet_prefix_only = true`
- `physical_chart_finite_bound_verified = false`
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

The actual hierarchy provider still has to answer every positive order demanded
by arbitrary common physical requests `J(M,P)` with one real assembled
coefficient object carrying theorem-backed derivative majorants, common compact
support, and exact retained-recurrence evidence under one coefficient state.
The current first-Picard/third-eta/A0-A1 value-row frontier (#319/#322/#326/#332)
does not yet provide that arbitrary-order assembled coefficient API.
