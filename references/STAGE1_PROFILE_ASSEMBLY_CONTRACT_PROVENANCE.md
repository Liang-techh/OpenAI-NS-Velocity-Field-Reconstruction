# Stage-1 profile-assembly identity contract provenance

## Scope

This increment is deliberately downstream of the coefficient-space work owned by
the dedicated Issue #1 backend lane. It does **not** implement
`AxisCoefficientSpace`, `naturalRemainder`, or a Picard iterate.

The integration baseline is repository `main`
`0ea9406195545d04cafd72fbed2ecf78e9f17ca8`, immediately after Agent 6 PR
#300 merged its fully-green angular actual-schedule eta-jet increment. That
upstream work remains partial: axial all-order jets, the certified infinite
weighted `AxisCoefficientSpace.ofJetFamily` bound, a complete typed reference /
fixed-point state, and genuine `naturalRemainder` are still absent. This
increment consumes none of those missing objects and does not duplicate Agent
6's backend lane.

The landed Stage-1 chain already owns the actual `SchedulePressure` datum, the
analytic low-|Z| cutoff choice, a common coefficient radius `epsilon`, the
wide theorem-selected `Lambda`, symbolic `C = exp(log_C)`, the reference pair,
and the scalar contraction certificate. A future fixed-point backend must not
silently mix a state from another admissible schedule or another theorem scale
with `NaturalProfileAssembly`.

`stage1_profile_assembly_contract.py` therefore creates one fail-closed identity
handshake from those existing owners. The handoff pins exactly:

- the outgoing schedule tuple `(P,m,lam,wait,h)` and `j`;
- the actual-schedule `sigma` and coefficient radius `epsilon`;
- wide `Decimal` `Lambda` and symbolic `log C`;
- the official coefficient window `[-11/10,11/10]`; and
- the natural amplitude convention `exp(Lambda*realPhase)/C`.

The contract uses exact equality, not tolerances or sampled agreement.

## Source mapping

The structure is tied to the official Lean commit pinned by
`references/provenance_manifest.json`, especially:

- `NavierStokes/AxisCoefficientSpace.lean`;
- `NavierStokes/AxisContraction.lean`;
- `NavierStokes/NaturalProfile.lean`; and
- `NavierStokes/SchedulePressure.lean`.

`actual_schedule_reference_axis_state` remains authoritative for
`sigma/epsilon`, while `diagnose_actual_schedule_scale_chain_wide` and
`NaturalPicardContractionCertificate` remain authoritative for `Lambda/log C`
and the scalar Picard gate. This increment does not duplicate their formulas.

## Verification

`tests/test_stage1_profile_assembly_contract.py` reconstructs the identity from
the real repository schedule fixture already used by the Stage-1 contraction
tests and independently compares it to the current reference-state and wide
scale providers. Regressions then perturb schedule identity, `j`, `sigma`,
`epsilon`, `Lambda`, `log C`, and the coefficient window one at a time and
require the handoff to fail closed. An alternate amplitude convention and a
non-`Decimal` wide scale are also rejected.

The production and regression blobs are unchanged from PR #302, whose Actions
run #2157 completed successfully across both full Python/NumPy jobs,
wheel/outside-checkout checks, and all slices. This replay still requires its own
fresh Actions on the exact post-#300 current main before merge.

These tests certify identity plumbing only. They do not supply coefficient
values, fit a profile, or constitute evidence for the fixed-point theorem.

## Truth boundary

Status remains `formal-structure`.

`full_reconstruction=false` and `paper_exact_velocity_available=false` remain
mandatory. Agent 6 still must complete the compatible coefficient backend and
`naturalRemainder`, the Picard/fixed-point state, and materialized
`phi/u/average/pressure`. This lane must then instantiate actual
`NaturalProfileAssembly`, close regular-inner-core / heat-exterior matching, and
establish support/moment/matching/cone certificates before Issue #1 can close.