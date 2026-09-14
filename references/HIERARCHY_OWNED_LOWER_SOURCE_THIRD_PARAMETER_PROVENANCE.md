# Section 5 hierarchy-owned third-eta `actualLowerSource` provenance

Status: **Stage 2 / formal-structure only**.

## Scope

This increment closes the next strict-lower derivative seam after the landed
hierarchy-owned third-eta angular preceding diffusion, axial preceding
diffusion, and regular Eq. (5.6) `Omega/X` rows.  The implementation baseline is
repository `main` commit
`95ba525c8aadab0ae4b05610b2042c9a02fec234`.

The new bridge is
`hierarchy_owned_lower_source_third_parameter_jet(...)` in
`background_repaired_history_source_third_parameter.py`.  It accepts only
`Section5LowerHistorySixthMixedHierarchy` and returns value through
`partial_eta^3 actualLowerSource` for one positive recursive order.

## Paper / implementation map

The strict-lower source is the source appearing in the PositiveAxis form of
Eqs. (5.3)-(5.7).  For recursive order `n`, its convolution rows use only
coefficient orders strictly below `n`; its two preceding-diffusion rows use
`phi_(n-1)` and `U_(n-1)`; and its regular vorticity row is
`Omega_(n-1)/X` from Eq. (5.6).

The new bridge composes already-landed hierarchy-owned data:

- repaired fifth-mixed `phi_j` jets;
- repaired sixth-mixed `U_j` jets and their fifth-mixed projections;
- fifth-mixed `beta_j=V_j/X` derived only through analytic Eq. (5.2);
- the analytic hierarchy-owned angular and axial third-eta
  `Z_(b-D)(Z_b F)` rows; and
- the analytic hierarchy-owned third-eta Eq. (5.6) `Omega/X` row.

The existing second-eta lower-source bridge remains authoritative for value,
first eta derivative, and second eta derivative.  Only the third derivative is
assembled here.

## Exact differentiation

Each quadratic strict-lower product is differentiated by the ordinary exact
third-order Leibniz rule

`(ab)''' = a''' b + 3 a'' b' + 3 a' b'' + a b'''`.

The displayed first-order operator `Z_b F` is differentiated analytically
through eta order three from the hierarchy-owned fifth-mixed profile jet.  The
two preceding-diffusion contributions and the `Omega/X` contribution are not
recomputed or supplied by callers: the bridge calls their landed hierarchy-owned
third-eta adapters.

Production therefore contains no finite-difference differentiation, fitted
coefficient family, generic cutoff substitution, sampled `V/X` division, or
caller-maintained `SourceJet` derivative table.

## Verification

`tests/test_background_repaired_history_source_third_parameter.py` constructs a
contiguous `Section5LowerHistorySixthMixedHierarchy` from coherent analytic
strong-jet providers.  The fixture is deliberately labelled as a regression
fixture rather than paper coefficient data.  The regression verifies:

1. exact projection of value/first/second rows to the landed second-eta bridge;
2. all four third-eta source components — angular, axial, pressure product and
   `Omega/X` — against centered eta differentiation of the older analytic
   second-eta path;
3. finite regular evaluation at `X=0`; and
4. fail-closed rejection of a non-strong hierarchy.

The centered finite difference is test-only and is not used in production or
promoted to proof evidence.

The compact Lemma 5.2 repair and the sixth-mixed-U/fifth-mixed-beta ownership
used by this composition are independently regression-tested in their landed
modules.  This increment does not replace those repairs with the analytic test
fixture.

## Truth boundary and workstream split

Genuine Issue-#1 order-zero strong profile data, positive-order unrepaired base
jets, moment/patch eta jets, and normalization inputs remain upstream contracts.
This increment closes only the hierarchy-owned `partial_eta^3 actualLowerSource`
composition seam.  It does **not** establish hierarchy-owned
`partial_eta^3 f_n`, `partial_eta^3 W_n^(0)`, `partial_eta^2 W_n^(1)`, a
materialized positive-order coefficient, Picard convergence, an all-order
coefficient family, paper-exact velocity, or full reconstruction.

Uniform/stagewise `C[j,m]` bounds, recursive SlowBorel/DiagonalScale cutoff
selection, arbitrary-target tail order, and all-order residual convergence are
kept outside this increment and remain in Agent 7's separate convergence lane.
No finite-prefix convergence result is promoted here.

Mandatory truth flags remain:

- `full_reconstruction=false`;
- `paper_exact_velocity_available=false`;
- Stage 2 status `formal-structure`.
