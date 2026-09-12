# Hierarchy-owned lower-source parameter provenance

## Scope

This increment combines the already-landed Section 5 derivative layers into one
analytic eta jet of the complete strict-lower `actualLowerSource`.  The
implementation is
`src/openai_ns_reconstruction/background_repaired_history_source_parameter.py`.

It consumes only a `Section5LowerHistoryPhiThirdMixedHierarchy`.  Ordinary
source values remain delegated to the previously landed
`positive_axis_point_data_from_lower_history` path.  The eta derivative is
assembled from the same hierarchy-owned coefficient data: Eq. (5.2) regular
flux jets, analytic derivatives of the angular/axial convolution terms,
analytic preceding-diffusion parameter jets, the pressure-product derivative,
and the hierarchy-owned analytic Eq. (5.6) `partial_eta(Omega/X)` bridge.
No production finite difference, sampled derivative table, generic cutoff, or
fitted coefficient is introduced.

## Paper connection

The relevant paper layer is Section 5, especially Eqs. (5.3)-(5.7), Lemma 5.2,
and the compact moment correction around Eq. (5.14) / Appendix A.  This is the
source-derivative layer needed before the already-landed identity
`partial_eta(G f_n) = G(partial_eta f_n)` can be used in a genuine next Picard
application.

For the strict-lower source at coefficient order n, the implementation
differentiates exactly the same rows as the landed value path:

- angular lower-order products and angular preceding diffusion;
- axial lower-order products and axial preceding diffusion;
- the lower pressure product; and
- `Omega_(n-1)/X` from Eq. (5.6).

The stronger phi/U/regular-flux derivatives are not supplied independently to
this bridge; they are queried through the strong hierarchy that already owns the
Lemmas-5.2 repaired derivative layers.

## Fail-closed boundary

A complete analytic source eta jet needs strong order-zero data.  In particular,
the preceding-diffusion row requires the leading phi third-mixed jet and the
Eq. (5.6) derivative requires leading fourth-mixed U data.  The real leading
profile remains owned by Issue #1.  If these data are absent, the hierarchy
fails closed rather than manufacturing them.

The following remain upstream or incomplete:

- the genuine Issue #1 leading mixed profile;
- unrepaired positive-order base mixed jets and moment/patch eta-jets;
- an analytic eta jet of the full Eq. (5.7) forcing `f_n` (including its
  explicit geometry factors after `actualLowerSource` is inserted);
- the genuine hierarchy-owned `k=1` Eq. (5.7) Picard application;
- recursive coefficient materialization, recursive cutoff-scale completion, and
  Proposition 5.3 all-jets truncation/residual decay;
- paper-exact background or final velocity.

Accordingly Stage 2 remains `formal-structure` and
`paper_exact_velocity_available=false`.

## Validation

`tests/test_background_repaired_history_source_parameter.py` constructs a
strong analytic order-zero fixture plus an actual positive-order
`Lemma52MomentRepair` compact five-bump repair.  At recursive order two it
compares every analytic source derivative (`angular`, `axial`,
`pressure_product`, and `omega_quotient`) against a centered finite difference
of the independently landed source-value path.  The regression is run both on
the axis (`X=0`) and inside an active compact-repair region.  The source value
itself must exactly equal the old path.

A separate regression removes leading fourth-mixed U data and verifies that the
new bridge fails closed instead of pretending the Eq. (5.6) eta derivative is
owned.  These are floating-point implementation checks, not proofs of the
paper's all-order estimates, uniform moment invertibility, recursive cutoff
schedule, or residual decay.
