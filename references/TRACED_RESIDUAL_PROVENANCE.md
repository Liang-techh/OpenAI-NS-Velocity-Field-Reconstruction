# CandidateFromLimits traced-residual bridge provenance

Status: **formal-structure only**. This note does not promote Stage 7 or make `paper_exact_velocity_available` true.

## Pinned source mapping

The repository pins OpenAI `NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The relevant official definitions are:

- `NavierStokes/SpacetimeEndpoint.lean`:
  `extendTrace T f L z = if z.1 < T then f z else L z.2`.
- `NavierStokes/CandidateFromLimits.lean`:
  `tracedResidual u p L` applies that `extendTrace` at `T=1` to
  `PastExtension.pastResidual u p`, using the degree-zero endpoint tensor
  `(L x 0).curry0` on the terminal/future branch.
- `NavierStokes/SpacetimeGluing.lean` and
  `NavierStokes/CandidateFromLimits.lean` then build the future side from the
  actual normal endpoint jets by Taylor--Borel extension.

## Executable artifact

`src/openai_ns_reconstruction/traced_residual.py` adds two deliberately narrow
helpers:

1. `traced_residual_from_full_spacetime_jets(...)` implements the exact
   value-level `extendTrace` branch. For `t<T` it evaluates the supplied
   closed-past residual; for `t>=T` it evaluates only the validated degree-zero
   full spacetime endpoint tensor. Branches short-circuit, so missing past data
   are not sampled on the terminal/future side and endpoint jets are not
   demanded on the open past.
2. `formal_force_from_full_spacetime_jets(...)` composes that trace with the
   already-landed `endpoint_jets.py` normal-time contraction and
   `endpoint_borel.py` locally finite right extension. Thus the executable
   branch structure matches the pinned construction at the level of values.

`tests/test_traced_residual.py` independently checks the `t<T` / `t>=T`
branch boundary, exact endpoint value, a finite Taylor--Borel future fixture,
future zero support from the cutoff scale, and fail-closed malformed/nonfinite
endpoint tensors. The fixture is structural test data, not a surrogate OpenAI
velocity or force.

## What this does **not** prove

The official `CandidateFromLimits` theorem assumes, for every derivative order,
locally uniform `t -> 1-` limits of the **actual** Navier--Stokes residual's full
spacetime derivatives. Merely supplying arrays through the executable API does
not establish those limits. Likewise, the supplied Borel scales are not the
missing analytic compact-template derivative majorants.

Therefore this increment does **not** establish:

- the actual residual full-jet family of the unresolved final local field;
- locally uniform convergence of those jets to the supplied endpoint tensors;
- closed-past or cross-endpoint `C^infinity` smoothness;
- the genuine SpatialBorel derivative-bound schedule for the actual jets;
- an independently verified final forcing;
- the paper-exact final velocity/forcing construction.

In particular, this bridge must not be used as a disguised `f=R; R-f=0`
verification. Stage 7 remains `formal-structure`, and the existing paper-exact
gate remains fail-closed until the real upstream field and endpoint analytic
certificates are supplied.
