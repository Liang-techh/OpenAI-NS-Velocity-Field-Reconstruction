# Finite cofinal physical-tail ladder provenance

## Scope

Issue #2, Agent 7 downstream/parallel all-order analytic-closure lane.

This increment is stacked directly on PR #354.  It does not implement Agent 2
hierarchy jets, Eq. (5.7) coefficient recurrence, or coefficient extraction.
It advances only the quantifier structure of the existing provider-owned
SlowBorel/DiagonalScale tail certificates.

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
  - `NavierStokes.SlowBorelBase.exists_physical_uncut_tail`

The pinned physical-tail theorem is quantified by one finite derivative order
and one requested real power, while the admissible-scale construction supplies
one recursive schedule.  Therefore an executable route toward all-jets
flatness needs both a cofinal target family and compatibility of its finite
prefixes with one schedule.

## Canonical cofinal target family

`background_physical_tail_cofinal_ladder.py` uses the exact family

`level k = (max physical derivative order k, target physical q-power k)`.

For an exact rational request `(M,P)`, the least canonical level dominating it
is computed without floating arithmetic as

`k = max(M, ceil(P), 0)`.

This is the only cofinality statement made by the increment.  It is arithmetic
on target indices; it is not a proof that the hierarchy provider exists at
arbitrary level.

## Finite materialization

`certify_finite_cofinal_physical_tail_ladder(...)` materializes levels `0..K`
from one `TargetDrivenExactMajorantProvider`.  Each level is built once through
the landed physical-tail prefactor chain.  Every adjacent pair is then checked
with PR #354's `FiniteIncreasingPhysicalTailCertificate`, preserving all of its
premises unchanged:

- literal shared exact `C[j,m]` majorants;
- literal shared common-support witnesses;
- literal shared exact retained recurrence identities;
- the same own-order coefficient artifact/provider links;
- the same executable normalized jet-bound rows;
- the same local cutoff-scale prefix; and
- the same recursive doubling-envelope prefix.

A stateful provider that changes a retained `C[1,0]` only on the third/deeper
request is rejected by regression.  This specifically checks extension
stability beyond a single pair.

A request whose canonical level exceeds `K` fails closed.  The finite ladder is
never extrapolated into an unavailable hierarchy order.

## Truth boundary

Verified here:

- exact rational canonical-target cofinality arithmetic;
- one finite canonical ladder `0..K`;
- literal shared schedule/evidence prefix coherence at every adjacent level;
- exact physical-power dominance for every derivative row in each materialized
  canonical level.

Not verified here:

- provider totality for all positive orders;
- theorem-owned replacement of Agent 2 PR #353's conditional `C_n`, radial
  extent, `rho`, or `rho_prime` inputs;
- one completed infinite coefficient hierarchy or diagonal schedule;
- a numerical physical-chart finite-bound constant;
- an actual PDE residual inequality;
- Proposition 5.3;
- all-jets flatness;
- super-algebraic residual convergence;
- paper-exact velocity; or
- full reconstruction.

A finite prefix of a cofinal target family is still finite evidence.  Passing
regression or CI cannot promote it to an infinite convergence theorem.

## Upstream requirement returned to Agent 2

The exact remaining coefficient-level handoff is a **total stable
assembled/repaired positive-order provider**.  For every order requested by the
canonical ladder, it must theorem-own the actual coefficient object/state,
`C[j,m]` derivative majorants, common compact support, exact retained recurrence
decomposition, and the analytic data replacing PR #353's conditional Eq. (5.8)
tail hypotheses.  Extending to a deeper level may not alter any previously
emitted object or evidence row.
