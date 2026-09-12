# Hierarchy-owned positive-axis forcing parameter provenance

## Scope

This increment closes the analytic eta-derivative layer between the already
landed hierarchy-owned strict-lower `actualLowerSource` and the exact Section 5
Eq. (5.7) positive-axis forcing `f_n`.  The implementation is
`src/openai_ns_reconstruction/background_repaired_history_forcing_parameter.py`.

The value path is delegated unchanged to the existing
`background_positive_axis.positive_axis_forcing`.  The source derivative is
obtained only from `hierarchy_owned_lower_source_parameter_jet`, which in turn
requires the strong repaired phi/U history.  The only new differentiation is
of the explicit `eta / ell(h,eta)` geometry in the sixth forcing row.  No
production finite difference, sampled derivative table, generic cutoff, or
finite-order fitted coefficient is introduced.

## Paper connection

The relevant layer is Section 5, Eqs. (5.3)-(5.7), after the recursive
strict-lower source has been assembled.  For

`pressureSource = C^-2 pressureProduct - omegaQuotient/2`,

its eta derivative is obtained by the same linear combination of the
hierarchy-owned source eta jet.  The first three forcing rows remain zero; the
remaining rows are differentiated analytically:

- `2 xi pressureSource`;
- `2 angular`; and
- `2 axial - 4 eta X pressureSource / ell`.

In particular the final row includes both `partial_eta pressureSource` and the
explicit derivative of `eta/ell`, with `ell = 1 - 2 h eta^2`.

This produces a genuine hierarchy-owned analytic `partial_eta f_n` whenever
the hierarchy owns all required strict-lower strong jets.  It is the missing
input needed by the already-landed identity
`partial_eta(G f_n) = G(partial_eta f_n)`.

## Fail-closed boundary

The interface accepts only a `Section5LowerHistoryPhiThirdMixedHierarchy` and
does not accept a caller-supplied source derivative.  The underlying source
bridge requires strong order-zero phi/U derivative data.  Until Issue #1
materializes the genuine leading mixed profile, production requests that lack
those data fail closed instead of fabricating them.

The following remain open:

- the genuine Issue #1 leading mixed profile;
- unrepaired positive-order base mixed jets and moment/patch eta-jets;
- wiring this hierarchy-owned `partial_eta f_n` into the singular inverse to
  materialize a hierarchy-owned `(G f_n, partial_eta G f_n)`;
- the next nontrivial Eq. (5.7) Picard application containing the landed exact
  `A0 W + A1 partial_eta W` terms;
- recursive coefficient materialization, recursive cutoff-scale completion,
  Proposition 5.3 all-jets truncation/residual decay, and paper-exact velocity.

Accordingly Stage 2 remains `formal-structure` and
`paper_exact_velocity_available=false`.

## Validation

`tests/test_background_repaired_history_source_parameter.py` reuses the strong
order-zero analytic fixture plus the actual positive-order
`Lemma52MomentRepair` compact five-bump repair.  At recursive order two, both on
the axis and inside active compact-repair support, it compares the new analytic
forcing eta jet against a centered finite difference of the independently
landed forcing-value path.  The forcing value itself must be bitwise identical
to the existing value implementation.

A separate regression removes the leading fourth-mixed U provider and verifies
that the forcing eta bridge inherits the lower-source fail-closed behavior.
Finite differences occur only in the test oracle; they are not part of the
production implementation.  These tests are floating-point implementation
checks, not proofs of all-order estimates, recursive cutoff admissibility, or
truncation/residual decay.
