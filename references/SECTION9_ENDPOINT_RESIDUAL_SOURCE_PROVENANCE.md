# Section 9 endpoint residual-source binding provenance

Status: **formal-structure**. This increment closes one provenance/identity hole on the Issue #4 endpoint path. It does not derive a new residual estimate, construct the actual Section 9 correction sequence, identify the true endpoint limit, or prove the infinite Borel extension smooth.

## Why this layer exists

The landed endpoint pipeline already separates several mathematically different claims correctly:

- `section9_endpoint_ladder_bridge.py` admits a finite ladder of independently justified uniform physical-time derivative witnesses and converts them to Section 10 endpoint majorants;
- `time_localization.py` checks that the whole majorant validity interval lies in the official `t >= 3/4` unit plateau of the Section 10 time switch, so no transition-collar correction terms are silently dropped;
- `section10_endpoint_borel_admission.py` evaluates dense pre-endpoint residual jets, encloses candidate endpoint tensors with the majorant tail budgets, freezes the candidate, and only then invokes the finite Borel-prefix right-jet algebra.

Before PR #140, however, the endpoint-majorant ladder and the `past_full_jets` evaluator entered the Section 10 candidate admission as independent arguments. A future caller could therefore attach majorants justified for residual source/revision A to dense residual jets evaluated from source/revision B. Each component could be internally valid while the combined endpoint enclosure was semantically meaningless.

PR #140 introduced `src/openai_ns_reconstruction/section9_endpoint_residual_source.py` to close that source-identity gap. This follow-up hardens the same gate: source identity now binds not only each ordered provenance string but also the corresponding ordered **evidence kind** carried by the admitted Section 9 majorant ladder. Reusing an identical provenance label while silently changing a bound from, for example, `paper-derived` to `certified-numerical` is rejected.

## Interface and checks

`Section9ResidualEndpointSourceWitness` records:

- one Section 9 stage and spatial window;
- a stable `(source_id, source_revision)` pair;
- theorem-only source evidence (`analytic-theorem` or `formal-theorem`) and nonempty source provenance;
- the exact ordered majorant-evidence kinds for endpoint degrees `0..N`;
- the exact ordered majorant-evidence provenance strings for endpoint degrees `0..N`;
- the dense pre-endpoint full-spacetime residual-jet evaluator;
- separate theorem assertions that the residual provider contract is the intended one, that the admitted majorants apply to this source revision, and that the dense past-jet evaluator is from the same source revision.

`Section9ResidualEndpointSourceBinding` then requires exact agreement with an already-admitted `Section9EndpointLadderBridgeRecord` in stage, spatial window, derivative count, every ordered evidence kind, and every ordered majorant provenance string. The existing finite ladder must remain formally ready. Sampled/fitted source evidence is rejected and no source identity is inferred from numerical values.

The source-bound candidate entry point intentionally has **no separate `past_full_jets` argument**. Once a binding is formed, the Section 10 trace/Borel-prefix path can only consume the evaluator frozen into that source/revision witness. This prevents accidental provider swapping after the majorant ladder has been bound.

## Relation to the Section 10 time-localization target

This increment does not add another arbitrary time window. It reuses the already-landed official Section 10 time switch and the existing `t >= 3/4` endpoint-localization transfer gate. The fixed spatial cutoff / analytic cutoff-gradient cross-validation is also left unchanged because that path is already covered by `spatial_localization.py` and its independent product-rule regression.

The purpose here is narrower: make the next genuine bridge from the actual Section 9 residual estimates into the `t=1` endpoint extension fail closed unless the majorant evidence metadata and the residual jets are demonstrably about the **same residual revision**.

## Independent regression

`tests/test_section9_endpoint_residual_source.py` uses analytic/formal fixtures only. It checks that:

- one stable source revision survives through the existing finite endpoint-candidate and Borel-prefix admission;
- ordered majorant provenance from a different residual revision is rejected;
- changing an ordered majorant evidence kind while keeping the same provenance text is rejected;
- blank per-degree evidence kinds are rejected before binding;
- stage or spatial-window mismatches are rejected;
- sampled/numerical source evidence, blank revisions, or a missing same-source theorem assertion are rejected;
- the source-bound production entry point exposes no parameter that can swap in an unrelated past residual-jet provider;
- successful admission leaves all actual-residual, endpoint-limit, infinite-Borel, smooth-forcing, and paper-exact flags false.

No fixture is claimed to be a sample of the manuscript blow-up solution.

## Remaining boundary

Still required before Issue #4 can cross the true `t=1` boundary:

- materialize the genuine Section 9 correction sequence / Eq. (9.21) residual source on the required physical late-time compact;
- derive machine-linked uniform estimates for `partial_t D^n R` from that exact source, to all orders needed by the endpoint theorem rather than finite caller-supplied fixtures;
- prove locally uniform endpoint limits and their uniqueness for the actual residual jets;
- provide the analytic template bounds / diagonal-scale argument needed for the infinite Borel right extension and all-order smoothness;
- only then reconstruct the compact smooth forcing and independently certify closure/smoothness, finite energy, divergence-free structure, compact support, and the blow-up path.

Therefore `source_majorants_derived_from_actual_residual_verified=false`, `actual_section9_sequence_verified=false`, `endpoint_limits_constructed=false`, and `paper_exact_velocity_available=false` remain mandatory.
