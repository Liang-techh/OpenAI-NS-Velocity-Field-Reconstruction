# Section 10 endpoint-candidate to finite-Borel admission provenance

Status: **formal-structure only**. This increment does not identify an actual Section 9 endpoint limit, does not prove the infinite Borel extension smooth, and does not change `paper_exact_velocity_available=false`.

## Pinned source

Official Lean repository: `openai/NavierStokesAndEuler`

Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant statements and interfaces:

- `NavierStokes/CandidateFromLimits.lean`: endpoint candidates are built from locally uniform limits of the full spacetime derivatives of the actual closed-past residual.
- `NavierStokes/SpatialBorelExtension.lean`: `rightExtension_right_jets` matches all prescribed normal endpoint jets for the infinite jointly smooth Taylor--Borel future extension.
- `NavierStokes/SpacetimeGluing.lean`: `smoothExtension_contDiff` glues the closed-past residual to that future extension through `t=1`.

The repository has not yet established the locally uniform limits required by `CandidateFromLimits` for the actual Section 9 residual. Therefore a caller-supplied endpoint tensor family cannot be treated as the true limit merely because it has the correct dense shape or because a future Taylor prefix can interpolate it.

## Executable increment

`src/openai_ns_reconstruction/section10_endpoint_borel_admission.py` connects two already-landed finite conditional interfaces without strengthening either hypothesis.

First, `section10_endpoint_trace_enclosures` uses a certified finite `Section10EndpointMajorantLadder` and a caller-supplied pre-endpoint residual-jet family to produce, for each degree `0..N`, a closed dense-tensor ball centered at `D^n R(t,x)`. Its radius is the integrated endpoint-tail budget from the independently supplied majorant.

Second, the new admission gate evaluates the proposed full endpoint tensor once for every degree, freezes those tensors, and requires each frozen tensor to lie inside its corresponding trace enclosure. If any degree fails, the function exits before the Borel scale is queried. Only an admitted frozen family is passed to `certify_section10_borel_prefix_right_jets`, so the left consistency check and right finite-prefix certificate cannot silently consume different candidate values.

The check deliberately uses `distance <= radius` with no numerical tolerance that could enlarge the analytic majorant budget.

## Independent regression

`tests/test_section10_endpoint_borel_admission.py` uses an analytic dense-jet family with a known endpoint. Its pre-endpoint fixture differs from the true endpoint by exactly `0.1(1-t)` in one dense coordinate, while an independently supplied `alpha=0`, `C=0.2` majorant gives endpoint radius `0.2(1-t)`. At `t=0.9` the test therefore checks the closed-form distance `0.01` against the independently computed radius `0.02` before checking the future prefix normal jets.

A rejection regression perturbs one degree beyond its trace ball and supplies a Borel-scale callback that raises immediately if called. The expected trace rejection occurs without querying that callback, demonstrating fail-closed ordering. Malformed dense tensors also fail before right-prefix construction.

## Boundary retained

Passing this gate is only a **necessary finite-order consistency condition** for a proposed endpoint family. It does not establish:

- that the supplied majorant ladder was derived from the actual Section 9 residual;
- that the admitted candidate is the actual locally uniform endpoint limit, or that it is unique;
- the analytic `templateBound` estimates and summable derivative tails required by `SpatialBorelExtension`;
- arbitrary-order or infinite `rightExtension_right_jets`;
- joint smoothness of the force through `t=1`;
- final compact smooth forcing, bounded kinetic energy, or blow-up closure.

Accordingly the aggregate certificate keeps `source_majorants_derived_from_actual_residual_verified=false`, `actual_section9_residual_limits_verified=false`, `endpoint_limit_uniqueness_verified=false`, `analytic_template_bounds_verified=false`, `infinite_borel_right_jets_verified=false`, `all_order_borel_smoothness_verified=false`, `smooth_compact_forcing_verified=false`, and `paper_exact_velocity_available=false`.