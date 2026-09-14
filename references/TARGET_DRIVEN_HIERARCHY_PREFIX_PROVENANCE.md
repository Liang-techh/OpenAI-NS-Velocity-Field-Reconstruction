# Target-driven hierarchy prefix provenance

Issue: #2 — Agent 7 downstream all-order analytic-closure lane.

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`

Revision: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant formal results:

- `NavierStokes.SlowBorelBase.exists_template_jet_bound`
- `NavierStokes.SlowBorelBase.exists_admissibleScales`
- `NavierStokes.SlowBorelBase.ordinary_tail_bound`
- `NavierStokes.DiagonalScale.exists_diagonal_scales`
- `NavierStokes.DiagonalScale.exists_uniform_tail_order`
- the exact retained-recurrence identities consumed by the existing Agent-7 cancellation layer.

The pinned `exists_admissibleScales` theorem derives, for every coefficient order `j` and derivative order `m`, a finite compactness bound `C[j,m]` from the actual smooth coefficient and then chooses one recursive doubling schedule. `ordinary_tail_bound` requires `m <= J+3` and yields the scalar ordinary-jet tail exponent `h*(J+1)-m`. `exists_uniform_tail_order` supplies the target-order logic behind choosing `J` large enough for an arbitrary prescribed q-power once the all-order hypotheses are genuinely available.

## Executable increment

`src/openai_ns_reconstruction/background_target_driven_hierarchy_prefix.py` removes one caller-supplied-prefix seam from the downstream lane.

For one requested ordinary derivative order `m` and exact target q-power `N`, it:

1. selects the minimal pinned truncation order `J`, also enforcing the theorem-side `m <= J+3` requirement;
2. queries `C[j,k]` only through one hierarchy-owned provider for every `1 <= j <= J`, `0 <= k <= j+2`;
3. feeds those owned rows into the already-landed recursive SlowBorel/DiagonalScale schedule;
4. queries exact retained recurrence identities only through the same provider for every `0 <= n <= J`;
5. requires one hierarchy id, source revision, and coefficient-state id across the bound and recurrence layers.

There is intentionally no API argument containing caller-supplied `C[j,m]` rows or a caller-supplied list of retained identities.

## Fail-closed frontier behavior

The current Agent-2 hierarchy still has a finite derivative/order frontier. If a requested target requires a coefficient order beyond that frontier, the provider failure propagates and no target certificate is returned. The implementation does not insert default-zero coefficients, hand-written bounds, sampled maxima, fitted coefficients, or synthetic recurrence cancellation in production.

The independent regression includes a finite provider frontier and verifies that a target requiring the next order fails at that exact missing coefficient row before recurrence certification begins. Separate regressions reject coefficient-state cross-wiring, target/derivative-loss cross-wiring, missing provider methods, and floating target powers.

## Truth boundary

This increment proves only a **finite target-driven prefix conditional on the provider successfully producing every requested hierarchy-owned artifact**. It does not prove that the provider is total for all orders, does not construct an infinite diagonal schedule, and does not certify the actual PDE residual tail.

Mandatory status remains:

- `finite_target_prefix_only = true`
- `all_order_hierarchy_verified = false`
- `infinite_diagonal_schedule_verified = false`
- `pde_residual_tail_verified = false`
- `all_jets_flat = false`
- `super_algebraic = false`
- `paper_exact = false`
- `full_reconstruction = false`
- `paper_exact_velocity_available = false`

The exact upstream request to Agent 2 is now executable rather than merely descriptive: expose the real hierarchy through the provider interface for every requested order, including hierarchy-owned derivative/support bound rows and exact retained recurrence identities under one source/coefficient state. Until that provider can answer arbitrary requested target orders, no infinite/all-order claim is permitted.
