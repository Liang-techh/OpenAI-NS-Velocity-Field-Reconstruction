# Section 9 derived-endpoint-majorant source-binding provenance

## Scope

This increment connects two already-landed formal-structure layers without adding a new numerical localization path:

1. `section9_endpoint_majorant_adapter.py`, which derives one bounded Section 10 endpoint majorant from a theorem-facing **uniform** Eq. (9.18) residual envelope on the official late-time support; and
2. `section9_endpoint_residual_source.py`, which binds a finite endpoint-majorant ladder and a pre-endpoint dense-jet provider to one stable residual-source revision.

`section9_endpoint_majorant_source_binding.py` requires a contiguous family of uniform residual-envelope witnesses for endpoint degrees `0..N`, checks that they share one source id/revision, Section 9 stage, and spatial window, derives every row through the existing analytic adapter, assembles the existing finite ladder record, and finally delegates to the existing residual-source binding gate. Evidence kind and provenance identity are therefore rechecked by the pre-existing source gate rather than reimplemented here.

## What this closes

Before this bridge, callers could run the uniform Eq. (9.18) adapter degree-by-degree and separately assemble a source-bound ladder. The missing connective layer made it easy for future integration code to accidentally mix revisions or bypass the intended ordering/provenance checks while wiring those two stages together.

The new adapter fails closed on:

- missing or duplicate endpoint degrees;
- mixed residual source revisions;
- mixed Section 9 stages or spatial windows;
- disagreement between the envelope ladder and the frozen source witness; and
- evidence-kind/provenance mismatch caught by the existing `Section9ResidualEndpointSourceBinding` gate.

Regression includes an explicit evidence-kind swap to verify that this bridge truly passes through the existing source gate rather than only checking source strings.

## Truth boundary

This remains `formal-structure` only. The input `Section9UniformResidualEnvelopeWitness` objects are still theorem-facing inputs. This increment does **not** derive the uniform Eq. (9.18) constants, derivative losses, logarithmic powers, or all-order flat-remainder bounds from the actual Eq. (9.21) residual construction. It does not prove the pre-endpoint dense jets are the actual manuscript residual except through caller-supplied theorem assertions already required by the source witness.

Accordingly the following remain false/unproved:

- `source_majorants_derived_from_actual_residual_verified`;
- `actual_section9_sequence_verified`;
- endpoint-limit construction/uniqueness;
- infinite Borel right-jet convergence and all-order smoothness;
- smooth compact forcing;
- bounded-energy and blow-up closure; and
- `paper_exact_velocity_available`.

The next substantive step is still to obtain machine-derived uniform Eq. (9.18) data from the genuine Eq. (9.21) residual revision and feed that data into this now source-locked path.
