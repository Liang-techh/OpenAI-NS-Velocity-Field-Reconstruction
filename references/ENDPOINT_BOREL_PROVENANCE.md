# Section 10 endpoint Taylor--Borel extension provenance

Status: **formal-structure**.  This file records an executable part of the
pinned endpoint-extension mechanism; it is not a claim that the final force has
been constructed.

Pinned source: `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

## Source mapping

- `NavierStokes/BorelExtension.lean`
  - `term b j v s = cutoff (b*s) * (s^j/j!) * v`;
  - the all-order extension is the sum of these terms;
  - the selected scale grows through `DiagonalScale.doublingEnvelope` and the
    Lean proof uses derivative-majorant bounds to establish C-infinity
    convergence and the prescribed jets at zero.
- `NavierStokes/DiagonalScale.lean`
  - `doublingEnvelope b 0 = max 1 (b 0)`;
  - `doublingEnvelope b (n+1) = max (b (n+1)) (2*doublingEnvelope b n)`;
  - doubling implies only finitely many cutoff supports reach any fixed
    positive distance from the endpoint.
- `NavierStokes/SpatialBorelExtension.lean`
  - `rightExtension T` evaluates the spatial Taylor--Borel extension at
    `(t-T,x)`;
  - it has the prescribed normal time jets at `t=T` and vanishes for
    `t >= T+1`.
- `NavierStokes/SpacetimeGluing.lean`
  - `glue T f g` uses the **closed-past branch at `t=T`** and the future branch
    only for `t>T`;
  - the final `smoothExtension` glues the closed-past field to that right
    extension after proving all normal jets match.
- `NavierStokes/CandidateFromLimits.lean`
  - the force construction additionally requires locally uniform left limits
    of every full spacetime derivative of the *actual* Navier--Stokes residual.

## Executable increment

`src/openai_ns_reconstruction/endpoint_borel.py` now provides:

1. the exact natural-number `doublingEnvelope` recurrence;
2. the cutoff-monomial Taylor coefficient `s^j/j!`;
3. a pointwise evaluator for the nominal infinite right-extension family.  For
   every finite floating `s=t-T>0`, evaluation stops at the first degree with
   `scale_j*s >= 1`; doubling guarantees all later terms are also outside the
   cutoff support;
4. endpoint value handling and the exact `t >= T+1` zero-support consequence;
5. an algebraic closed-past/future glue helper whose documentation explicitly
   does **not** assert smooth matching;
6. `close_left_open_past(...)`, an execution adapter for the common situation
   where an upstream residual evaluator is available only for `t<T` while its
   endpoint trace is supplied separately by the degree-zero endpoint jet.  It
   leaves `t<T` unchanged, fills `t=T` with `jet(0,x)`, and fails closed for
   `t>T`; `glue_to_left_open_past(...)` then reuses the existing closed-past
   branch convention rather than silently switching the join to the future
   branch.

Tests independently sum the active terms, check the doubling recurrence,
recover first/second right jets for a manufactured finite jet family, verify
future support without touching unavailable jets, and exercise fail-closed
input validation.  The left-open completion regression uses a past callable
that deliberately raises at `t>=T` and deliberately does **not** match the
supplied endpoint value; this confirms both that the adapter never samples the
unavailable past field at the join and that the adapter itself does not claim
continuity.  These fixtures test the adapter only and are not profile or force
surrogates.

## Boundary that remains open

The runtime cutoff is the repository's explicit C-infinity representative with
the same plateau/support geometry as the pinned `SmoothCutoffs.cutoff`; no
pointwise equality with Mathlib's noncomputable transition-collar values is
claimed.

More importantly, the executable evaluator accepts `jet(j, x)` and
`local_scale(j)` as inputs.  `close_left_open_past(...)` only materializes a
closed-past *candidate value* at the endpoint; it is not evidence that the
incoming residual converges to that value.  The implementation does **not** yet
derive the spatial template bounds that justify the official all-order scale
choice, prove the derivative majorants, establish locally uniform convergence
of the actual residual jets as `t -> 1-`, or prove that all normal jets of the
closed-past completion match the right extension.  Those are precisely the
hypotheses needed before `SpacetimeGluing.smoothExtension` can be instantiated
for the real candidate.  The unresolved upstream local field from Issues
#1--#3 is also still required.

Therefore this increment must not change
`paper_exact_velocity_available=false`, must not upgrade Stage 7 beyond
`formal-structure`, and must not be used to manufacture a final forcing by
setting `f` equal to a numerically reconstructed residual.
