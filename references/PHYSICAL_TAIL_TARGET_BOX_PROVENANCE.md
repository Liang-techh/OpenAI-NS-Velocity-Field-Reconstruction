# Finite physical-tail target-box provenance

## Scope

Issue #2, Agent 7 downstream/parallel all-order analytic-closure lane.

This increment is stacked directly on PR #358. It does not implement Agent 2
hierarchy jets, Eq. (5.7) Picard recurrence, compact moment repair, or
coefficient extraction. It takes one finite quantifier step on top of the
already provider-owned, prefix-coherent physical-tail ladder.

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`

Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant source:

- `NavierStokes/DiagonalScale.lean`
  - `NavierStokes.DiagonalScale.doublingEnvelope`
  - `NavierStokes.DiagonalScale.exists_diagonal_scales`
  - `NavierStokes.DiagonalScale.exists_uniform_tail_order`
- `NavierStokes/SlowBorelBase.lean`
  - `NavierStokes.SlowBorelBase.exists_admissibleScales`
  - `NavierStokes.SlowBorelBase.physical_composition_bound`
  - `NavierStokes.SlowBorelBase.exists_physical_uncut_tail`

The formal all-order direction requires that every fixed physical derivative
order and every requested algebraic decay power can be reached by a sufficiently
late prefix of one recursively selected schedule. PR #358 records only a finite
prefix of a positive cofinal diagonal target family and explicitly does not
infer provider totality.

## Finite rectangular quantifier

`background_physical_tail_target_box.py` certifies one finite rectangle

`0 <= m <= M`, `1 <= N <= P`.

It selects no new coefficient data. Instead it requires PR #358's already
materialized cofinal ladder to reach the single common canonical frontier

`K = max(M, P, 1)`.

The frontier certificate is the literal level-`K` object already stored in that
ladder. Every target cell `(m,N)` is bound to the exact derivative-`m`
`PhysicalTailPrefactorRow` of this same frontier certificate. The cell also
records its least canonical level `max(m,N,1)`, but the proof obligation for the
whole box is discharged by the one common frontier prefix.

Consequently every certified cell inherits, without caller substitution:

- the same provider-owned exact `C[j,m]` prefix;
- the same common compact-support evidence;
- the same exact retained recurrence identities;
- the same own-order coefficient dependency links;
- the same recursive SlowBorel/DiagonalScale schedule prefix;
- the same first-omitted-order ledger and uncut radius;
- the same finite physical-chart bound shape; and
- the exact derivative-row prefactor `m! * 2^-J * D^m` together with a physical
  q-power at least the requested integer `N`.

The implementation rejects a box whose frontier exceeds the materialized
ladder. It also rejects cross-wired frontier objects rather than accepting an
equal-looking certificate supplied independently.

## Regression

The focused regression materializes a three-level canonical ladder and requests
`M=2`, `P=3`. One literal level-3 frontier must cover all nine cells in the
`3 x 3` finite rectangle. The regression checks the exact minimal level of every
cell, common frontier identity, derivative-row alignment, common dyadic factor,
`D^m` exponent, physical-power dominance, finite-frontier failure, malformed
requests, and the non-promotion truth flags.

The provider used in regression is analytic fixture data, not the paper's
hierarchy. No small numerical residual or finite plot is interpreted as an
identity or convergence theorem.

## Truth boundary

Verified here:

- one finite rectangular quantifier over physical derivative orders and positive
  integer algebraic powers;
- one common canonical frontier prefix for every cell in that rectangle;
- exact binding of each cell to the frontier derivative-row prefactor and
  physical exponent; and
- fail-closed behavior beyond the finite cofinal-ladder frontier.

Not verified here:

- provider totality for arbitrary positive coefficient order;
- theorem-owned replacement of Agent 2 PR #353's conditional `C_n`, radial
  extent, `rho`, or `rho_prime` inputs;
- one completed infinite coefficient hierarchy or diagonal schedule;
- numerical materialization of the physical-chart existential constant `D`;
- an actual PDE residual inequality;
- Proposition 5.3;
- all-jets flatness;
- super-algebraic residual convergence;
- paper-exact velocity; or
- full reconstruction.

A finite rectangle, even one dominated by a single coherent prefix, is not an
infinite all-jets/super-algebraic proof.

## Upstream requirement returned to Agent 2

The coefficient-level blocker remains a total stable assembled/repaired
positive-order provider. For every order demanded by increasingly large finite
target boxes, it must theorem-own the actual coefficient object/state,
`C[j,m]` derivative majorants, common compact support, exact retained recurrence
decomposition, and the analytic data replacing PR #353's conditional Eq. (5.8)
majorant inputs. Increasing the requested frontier must preserve every already
emitted coefficient object and evidence row literally.