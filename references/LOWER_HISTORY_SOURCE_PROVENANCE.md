# Section 5 lower-history PositiveAxis source provenance

Status: **formal-structure**.

## Pinned official source

This increment is mapped directly to
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/PositiveAxisSystem.lean`
  - `lowerConvolution`
  - `angularConvection`
  - `axialConvection`
  - `lowerSource`
  - `precedingDiffusion`
  - `actualLowerSource`
  - `baseAtOrderZero`
- `NavierStokes/PositiveAxisExistence.lean`
  - `lowerHistoryData`
- `NavierStokes/SlowRecursion.lean`
  - the induction bridge from already-constructed lower profiles to
    `PositiveAxisExistence.lowerHistoryData`.

For positive order `n`, the official source contains only orders strictly below
`n`.  The convective and pressure-product sums are therefore
`i=1,...,n-1`, `j=n-i`; the known preceding viscosity is

`Z_(b-D) Z_b F_(n-1)`, with `b=power+lambda_(n-1)`,

and the radial residual source is the smooth quotient `Omega_(n-1)/X`.

## Executable artifact

`src/openai_ns_reconstruction/background_lower_history_source.py` adds:

- `ProfileSecondJet`, containing exactly the pointwise second derivatives needed
  to evaluate the preceding double-`Z` viscosity analytically;
- `preceding_diffusion_from_second_jet`, an analytic evaluation of the pinned
  `Z_(b-D) Z_b` expression, with no numerical differentiation in production;
- `positive_axis_point_data_from_lower_history`, which constructs the order-zero
  `PositiveAxisBaseJet` and the full `PositiveAxisSourceJet` from only the
  already-constructed orders `< n`;
- automatic construction of the `Omega_(n-1)/X` source through the already
  landed Eq. (5.6) regular quotient evaluator, interpreting `beta_j=V_j/X`.

This removes the need to hand-supply the four `SourceJet` scalars when the lower
coefficient history is available pointwise through second jets.  It is the
executable counterpart of the official `actualLowerSource`/`lowerHistoryData`
packaging used before the positive-order Eq. (5.7) solve.

## Independent checks

`tests/test_background_lower_history_source.py` does not compare the new double
`Z` routine with itself.  It evaluates the inner `Z_b` directly on an explicit
polynomial and obtains the outer radial/parameter derivatives by centered finite
differences, then compares that independently composed `Z_(b-D) Z_b` value with
the production analytic formula.

A separate order-2 regression expands the unique strict-lower convolution term
`i=j=1` by hand for the angular, axial and pressure-product sources, and checks
the `Omega_1/X` wiring against the previously independently tested Eq. (5.6)
primitive.  Order 1 is also exercised with exactly order-zero history, proving
that no current-order coefficient is accessed by this adapter.

## Boundary / remaining blocker

This does **not** materialize the lower profiles.  The second jets supplied to
the adapter must ultimately come from the genuine recursively solved and
Lemma-5.2-repaired coefficient hierarchy; sampled/fitted jets are not a
paper-exact substitute.  The current Stage-1 work has only begun to materialize
reference-pair coefficient functions and has not yet produced the complete
Theorem 4.6 fixed point required for the first genuine Stage-2 base profile.

This increment also does not certify a common complex strip, the Lemma-5.1
Picard constants `C_n`, convergence of the actual positive-order Picard series,
uniform nonvanishing for the compact moment repair, the infinite cutoff
schedule, or Proposition 5.3 all-jets-flat residual decay.

Therefore Stage 2 remains **formal-structure** and
`paper_exact_velocity_available` remains `false`.
