# Section 10 finite Borel-prefix right-jet provenance

Status: **formal-structure only**. This increment does not make the velocity or force paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

Official Lean repository: `openai/NavierStokesAndEuler`

Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant statements:

- `NavierStokes/SpatialBorelExtension.lean`
  - `rightExtension_contDiff`
  - `rightExtension_time_jets`
  - `rightExtension_right_jets`
- `NavierStokes/SpacetimeGluing.lean`
  - `smoothExtension`
  - `smoothExtension_contDiff`

The pinned `rightExtension_right_jets` theorem states that the one-sided time jets of the full jointly smooth Taylor--Borel right extension at `t=T` equal the prescribed coefficient family. `SpacetimeGluing.smoothExtension_contDiff` uses exactly that result to match the future normal jets to the closed-past normal trace.

## Executable increment

`src/openai_ns_reconstruction/section10_borel_prefix_right_jets.py` implements only the finite-prefix algebra underlying that theorem.

For a requested finite degree `N`, it uses the same landed `BorelRightExtension` doubling-envelope scales `a_j`. For every `n<=N`, it records the exact common right-neighborhood radius

`rho_n = 1 / (2 max_{j<=n} a_j)`.

On `0 <= t-T <= rho_n`, every cutoff in the finite prefix `j=0,...,n` lies on the unit plateau because `a_j (t-T) <= 1/2`. Thus that prefix is exactly

`sum_{j=0}^n (t-T)^j/j! * a_j(x)`,

and its degree-`n` endpoint time derivative is exactly `a_n(x)`. The coefficient `a_n(x)` is obtained by contracting the validated dense full spacetime endpoint tensor in the time direction in every derivative slot; it is not accepted as a separate caller vector.

The plateau radius is stored as an exact `Fraction`, so a very small certified neighborhood is not rounded to a fake binary64 zero.

## Independent regression

`tests/test_section10_borel_prefix_right_jets.py` uses a finite-support analytic endpoint-jet family. Production computes the endpoint derivative by exact prefix algebra. The regression instead evaluates the already-landed Borel future branch at four positive times strictly inside the common plateau, solves an independent Vandermonde interpolation problem, and differentiates the recovered cubic at the endpoint. The recovered derivatives agree with the certified prescribed normal jets.

A separate test checks the exact common plateau radii for a nontrivial doubling schedule, while malformed dense tensor shapes and a non-Section-10 endpoint fail closed.

## Boundary retained

This increment deliberately does **not** claim the full Lean theorem `rightExtension_right_jets`. Passing it does not establish:

- that the supplied full endpoint tensors are locally-uniform limits of the actual Section 9 residual derivatives;
- the global analytic `templateBound` estimates used by `SpatialBorelExtension`;
- summable derivative tails or justification for interchanging an infinite Borel sum with endpoint differentiation;
- arbitrary-order / all-order right-jet matching;
- joint smoothness of the glued force through `t=1`;
- compact smooth forcing, bounded kinetic energy, or blow-up closure.

Those remain upstream/downstream proof obligations. In particular, `infinite_borel_right_jets_verified=false`, `all_order_borel_smoothness_verified=false`, `smooth_compact_forcing_verified=false`, and `paper_exact_velocity_available=false` are hard-coded truth boundaries of the returned aggregate certificate.
