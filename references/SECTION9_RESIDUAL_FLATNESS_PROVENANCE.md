# Section 9 residual-flatness theorem admission provenance

## Increment

This is a single fail-closed Section 9 closure increment for Issue #3.  It does
not build another surrogate correction stage.  Instead it pins the exact
formal theorem that converts sufficiently strong finite-stage improvement data
into all-order flatness of the limiting physical Navier--Stokes residual.

Pinned formal source:

- repository: `https://github.com/openai/NavierStokesAndEuler`;
- commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`;
- file: `NavierStokes/DiagonalResidual.lean`;
- conclusion theorem:
  `NavierStokes.DiagonalResidual.allJetsFlat_residual_of_stages`.

The pinned proof uses
`NavierStokes.DiagonalResidual.residual_jetRate_of_stages`, which in turn uses
`NavierStokes.DiagonalResidual.residualDifference_jetRate`.  The formal source
states the key diagonal principle explicitly: for every requested derivative
order and residual decay power, one may choose a sufficiently advanced finite
stage.  A single fixed tail is not assumed to be flat to all orders.

## Exact handoff contract

`src/openai_ns_reconstruction/section9_residual_flatness.py` admits only an
external `lean-formal-export` application that identifies the exact domain,
filter, scale `q`, limiting velocity/pressure, stage velocity/pressure
families, divergent gain schedule, and the three loss schedules.  The export
must separately certify all hypotheses consumed by the theorem:

- the open/eventual domain and eventual `0 < q <= 1` conditions;
- smooth limit and stage fields;
- `g -> +infinity`;
- background stage jet rates;
- velocity and pressure tail jet rates through the requested finite prefixes;
- finite-stage residual jet rates; and
- the final theorem application itself.

These objects remain opaque stable identities.  Python does not invent an
encoding for Lean filters or fields, does not fit a gain function, and does not
sample a residual to infer an exponent.

## Why this is useful for Section 9

This gives the repository a precise convergence target for the residual-
improvement iteration.  Future executable/formal Section 9 work can export its
actual finite-stage record against this interface.  Once the real stage family
and rate witnesses exist, the formal source already supplies the quantifier
structure needed to pass from stage-by-stage improvement to an all-jets-flat
limiting residual.

This increment intentionally does **not** claim that the paper's Section 9
stage family, Eq. (9.21), correction amplitudes, or gain/loss schedules have
been materialized or replayed in this repository.

## Truth boundary

The admission remains `formal-structure` only:

- `all_jets_flat_theorem_machine_replayed=false`;
- `section9_iteration_machine_materialized=false`;
- `stage_fields_materialized=false`;
- `stage_rate_bounds_materialized=false`;
- `limit_velocity_materialized=false`;
- `navier_stokes_residual_values_materialized=false`;
- `paper_exact_velocity_available=false`.

No sampled/fitted/numeric-scan evidence is accepted and no surrogate wave,
background, or residual-improvement constant is introduced.
