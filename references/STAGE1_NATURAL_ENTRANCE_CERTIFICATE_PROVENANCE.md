# Stage-1 NaturalEntrance scalar-premise certificate provenance

## Scope

This increment is downstream-only and starts from repository `main`
`3cc438f29524acb068e40bf3a99a8293c9222628`.

It deliberately does **not** implement `AxisCoefficientSpace`, coefficient
products, `naturalRemainder`, `x2`, a fixed point, or materialized
`phi/u/average/pressure`. Agent 6 retains ownership of that upstream backend
lane; open PR #323 currently materializes the genuine first-Picard
`a^2 * phi1^2` natural-source constituent and is not duplicated here.

The landed Stage-1 scalar chain already provides actual-schedule
`remainderBound/remainderLip`, wide `Lambda`, and the Picard contraction gate.
The next downstream NaturalEntrance proofs require the same remainder constant
`K` through two exact scalar premises:

1. `AxisReference.stabilityScale epsilon K <= Lambda`; and
2. `||x-referencePair|| <= K/(2*Lambda)` for the genuine coefficient-space
   fixed point.

`stage1_natural_entrance_certificate.py` binds those premises without creating
or approximating the missing state.

## Pinned formal sources

Official formal source revision:
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant definitions/theorems:

- `NavierStokes/AxisReference.lean`
  - `stabilityScale`;
  - `positive_of_uniformMixedError`;
  - `log_slope_of_uniformMixedError`;
- `NavierStokes/AxisContraction.lean`
  - `exists_unique_natural_fixedPoint`;
  - `uniform_natural_fixedPoint`;
- `NavierStokes/NaturalEntrance.lean`
  - `coefficient_phi_lower`;
- `NavierStokes/NaturalProfile.lean`
  - `profileErrorConstant`.

For the pinned `R=5`, parameter-derivative order zero jet sums,
`stabilityScale` reduces exactly to
`1 + (14000/9) * remainderBound`. The existing wide scale selector already
implements this theorem-side algebra with upward-rounded 96-digit `Decimal`.
The fixed-point theorem supplies the norm-error radius
`remainderBound/(2*Lambda)`, exposed by the existing
`NaturalPicardContractionCertificate`.

## New certificate

`NaturalEntranceScaleCertificate` accepts only a
`WideNaturalScaleSelection` whose stored stability threshold is exactly the
pinned wide value and whose `Lambda` dominates it. The fixed-point error radius
is cross-bound to the landed Picard certificate; an eventual backend error
witness must be a finite nonnegative `Decimal` no larger than that radius.

The exact formal constants are recorded as rational `Fraction` values:

- entrance scaled radius `41/10`;
- strict coefficient-profile lower target `1/8`;
- log-slope evaluation radius `Y=4`;
- required `chi >= 99/100`;
- strict log-slope lower target `23/10`.

The log-slope chi gate rejects float/Decimal approximations. It accepts only
exact rational bounds meeting `99/100 <= chi <= 1`, so sampled `0.99` cannot
be promoted to the theorem hypothesis.

## Verification boundary

`tests/test_stage1_natural_entrance_certificate.py` includes an actual-schedule
regression using the same theorem-admissible `TailData(P=2,m=1,lam=0.05,
wait=30,h=0.01), j=0.05` fixture used by the landed Picard tests. It checks the
actual wide chain closes the NaturalEntrance scalar scale gate while
`paper_exact` and `full_reconstruction` remain false. Independent fail-closed
regressions reject a nonpinned stability threshold, `Lambda` below stability,
a backend norm error above `K/(2*Lambda)`, non-Decimal error metadata, and
inexact/insufficient chi bounds.

These tests do not prove that a coefficient-space fixed point exists in the
Python reconstruction; they only verify the exact scalar interface that the
future genuine backend must satisfy.

## Remaining boundary

Still required before any profile-level entrance claim can be promoted:

- Agent 6's genuine compatible coefficient state and complete
  `naturalRemainder` / `x2` / fixed-point materialization;
- a backend-owned norm-error witness tied to that exact state and reference
  pair;
- actual `phi/u/average/pressure` materialization;
- `NaturalProfileAssembly` instantiation;
- analytic/formal identification of the `chi >= 99/100` log-slope parameter
  region for the materialized profile;
- regular inner-core / heat-exterior matching and support/moment/matching/cone
  closure.

Status remains `formal-structure`. `full_reconstruction=false` and
`paper_exact_velocity_available=false` remain mandatory.
