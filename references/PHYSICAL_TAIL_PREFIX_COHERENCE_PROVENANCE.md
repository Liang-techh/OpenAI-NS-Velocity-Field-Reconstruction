# Physical-tail prefix-coherence provenance

## Scope

Issue #2, Agent 7 downstream/parallel all-order analytic-closure lane.

This increment checks one necessary finite coherence property of the paper's
all-order SlowBorel/DiagonalScale construction: a stronger physical
derivative/decay request must extend the *same* recursive schedule and the same
hierarchy-owned coefficient evidence already used by a weaker request.

It is stacked directly on PR #349.  It does not implement Agent 2's hierarchy
jets, Eq. (5.7) Picard recurrence, or coefficient extraction.

## Pinned formal source

Repository:

- `openai/NavierStokesAndEuler`
- commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant source:

- `NavierStokes/DiagonalScale.lean`
  - `NavierStokes.DiagonalScale.doublingEnvelope`
  - `NavierStokes.DiagonalScale.exists_diagonal_scales`
  - `NavierStokes.DiagonalScale.exists_uniform_tail_order`
- `NavierStokes/SlowBorelBase.lean`
  - `NavierStokes.SlowBorelBase.exists_admissibleScales`
  - `NavierStokes.SlowBorelBase.exists_physical_uncut_tail`

`exists_diagonal_scales` first chooses local positive scales and then applies
one recursive `doublingEnvelope`.  The final SlowBorel construction uses one
such schedule.  Request-dependent tail prefixes therefore must be compatible
with one schedule; independently rebuilt prefixes are not enough unless their
shared part is literally identical.

## Executable certificate

`background_physical_tail_prefix_coherence.py` builds two
`FinitePhysicalTailPrefactorCertificate` values from the **same**
`TargetDrivenExactMajorantProvider`.

The later request must dominate the earlier request in both finite derivative
budget and requested physical q-power.  The certificate then requires literal
prefix equality of all shared:

1. exact hierarchy-owned `C[j,m]` majorants;
2. common-support witnesses;
3. exact retained recurrence identities;
4. own-order artifact/provider dependency links;
5. executable normalized jet-bound rows;
6. local cutoff scales; and
7. recursive doubling-envelope scales.

It also checks the exact increasing-truncation arithmetic.  If the two selected
orders are `J0 <= J1`, every shared ordinary and physical derivative exponent
must gain exactly

`h * (J1 - J0)`,

the raw first-omitted fixed-stage power must gain exactly

`2 h * (J1 - J0)`,

and the dyadic tail prefactor ratio must be exactly

`2^(-(J1-J0))`.

The later uncut radius is required not to exceed the earlier radius.

## Fail-closed behavior

The regression includes a stateful provider that returns a different exact
`C[1,0]` majorant on the second request while keeping the same hierarchy,
source revision, and coefficient-state labels.  Each finite request would be
well-formed in isolation, but the pair is rejected because it cannot represent
two prefixes of one fixed recursive schedule.

No coefficient row, support table, recurrence identity, cutoff scale, omitted
order, or tail exponent is accepted directly from the caller.

## Truth boundary

This is a **two-request finite coherence certificate only**.

It does not establish:

- total provider availability for arbitrary positive order;
- one completed infinite hierarchy or infinite diagonal schedule;
- a numerical physical-chart constant;
- an actual PDE residual norm bound;
- Proposition 5.3;
- all-jets flatness;
- super-algebraic residual convergence;
- paper-exact velocity; or
- full reconstruction.

In particular, exact prefix coherence for two (or any finite number of)
requests is not an infinite convergence proof.

## Upstream requirement returned to Agent 2

The current Agent 2 frontier ends at a finite hierarchy-owned Eq. (5.7) Picard
derivative triangle.  The next coefficient-level handoff needed by this lane is
a total assembled/repaired positive-order coefficient provider.  For every
requested order it must return one stable actual coefficient object/state whose
theorem-backed `C[j,m]` majorants, common compact support, and exact retained
recurrence evidence remain identical when a later request extends the prefix.

A provider whose previously emitted coefficient or bound rows change when a
deeper order is requested cannot instantiate one paper SlowBorel schedule and
must fail closed.
