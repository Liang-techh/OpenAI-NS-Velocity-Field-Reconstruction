# Target-driven hierarchy support provenance

Issue: #2  
Lane: Agent 7 all-order analytic closure  
Parent: PR #315 (`background_target_driven_hierarchy_prefix.py`)

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`  
Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant formal facts:

- `NavierStokes.AssembledSlowBase.extendedCoefficient_support`
- `NavierStokes.AssembledSlowBase.extendedCoefficient_compactSupport`
- `NavierStokes.SlowBorelBase.exists_template_jet_bound`
- `NavierStokes.SlowBorelBase.exists_admissibleScales`
- `NavierStokes.DiagonalScale.exists_uniform_tail_order`

The pinned assembled coefficient theorem proves, for every positive coefficient
order and all four coefficient components, a common support enclosure

`tsupport(extendedCoefficient[n,i]) ⊆ [-1, B^2/2] × [-outer, outer]`,

where the same scheme `B` and the same `commonWindow.outer` are used across the
positive-order sequence.  Compact support follows from this enclosure.  The
SlowBorel template derivative bounds are then obtained from smoothness on a
fixed compact inner-coordinate set; they are not numerical assumptions on the
summed output.

## Executable increment

`background_target_driven_support_prefix.py` extends the target-driven provider
contract without accepting caller-supplied support rows.  For exactly the finite
order `J` selected by the pinned target-tail arithmetic, the provider must emit
one `HierarchyCoefficientSupport` per positive coefficient order.

Each support witness is fail-closed on all of the following:

1. one hierarchy id, source revision, and coefficient-state id;
2. exact `Fraction`/integer geometry only -- binary floats are rejected;
3. the pinned radial box `[-1, B^2/2]` and symmetric parameter box
   `[-outer, outer]`, with `B>0` and `outer>1`;
4. exactly the four assembled coefficient components `phi`, `axial`, `beta`,
   and `pressure`;
5. the pinned support and compact-support theorem names;
6. nonfuture hierarchy dependencies containing the coefficient's own order;
7. one literally common support box across every positive order in the selected
   prefix; and
8. at least one identical own-order `(artifact, provider)` dependency shared
   with **every** `C[j,m]` row at that coefficient order.

That final check prevents a support theorem for one coefficient object from
being silently paired with derivative bounds for a different object merely
because both carry the same hierarchy label.

## Deliberate boundary

This increment proves only a finite, target-selected support/derivative/
recurrence prefix.  It does **not** manufacture the missing support data for the
current finite Agent-2 hierarchy, and it does not prove that the provider can
answer arbitrary orders.  In particular:

- `finite_target_prefix_only = true`
- `all_order_hierarchy_verified = false`
- `infinite_diagonal_schedule_verified = false`
- `pde_residual_tail_verified = false`
- `all_jets_flat = false`
- `super_algebraic = false`
- `paper_exact = false`
- `full_reconstruction = false`
- `paper_exact_velocity_available = false`

No finite support list, finite plot, small residual, or numerical near-zero is
promoted to an infinite convergence theorem.

## Exact upstream handoff to Agent 2

The current real hierarchy reaches the landed `partial_eta^3(Omega_k/X)` layer
but still has a finite derivative/order frontier.  To discharge this new gate,
Agent 2's eventual arbitrary-order coefficient source must expose, for each
requested positive order, the actual assembled coefficient object together with
its `extendedCoefficient_support` / compact-support provenance under the same
hierarchy/source/coefficient-state used by its analytic `C[j,m]` rows and exact
retained recurrence decomposition.  The support witness must refer to the same
own-order artifact as those derivative-bound rows; an unrelated box or a
caller assertion is intentionally insufficient.
