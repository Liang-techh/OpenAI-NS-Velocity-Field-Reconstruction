# Section 9 finite-prefix evaluation provenance

Status: **formal-structure only**.

This increment does not construct the actual Section 7/8 correction fields, does not identify a runtime callable with the manuscript's fixed cutoff pointwise, does not construct the infinite Eq. (9.21) sums, and does not change `paper_exact_velocity_available=false`.

## Primary source

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026), Lemma 5.4 and Section 9 / Proposition 9.9, especially Eqs. (9.20)–(9.21). The canonical source is recorded in `references/SOURCES.md`.

The reconstruction plan records the Section 9 assembly schematically as

- `A = A_0 + sum_j chi(a_j q) A_j`,
- `B e_theta = B_0 e_theta + sum_j chi(a_j q) B_j`,

with the pressure correction summed analogously in Eq. (9.21). The scale sequence is theorem-selected; it must not be inferred from finite numerical fitting.

## Upstream facts reused

This module deliberately reuses, rather than reimplements, the landed structural gates:

1. `section9_correction_extension_admission.py` admits one positive stage only after complete `A_j/B_j/p_j` theorem/proof witnesses agree on a common `q_big`, an actual schedule scale `a_j`, the cutoff-support fact, and smooth zero extension.
2. `section9_local_sum_admission.py` records the independent infinite-doubling-schedule hypothesis and its local-finiteness arithmetic.
3. `section9_stage_physical_support.py` already bridges admitted scales to whole physical slabs. The present increment does not duplicate that support bridge.

## What this increment adds

`section9_finite_prefix.py` adds an evaluation-only adapter for a **contiguous finite positive-stage prefix** `j=1,...,N` at one exact rationalized `q` with `0<q<q_big`.

It requires:

- one complete `Section9CorrectionStageAdmissionCertificate` for every stage `1,...,N`,
- exact agreement with the supplied global common domain and first schedule scale,
- exactly one finite numerical `A_j` three-vector, `B_j` scalar, and `p_j` scalar for every admitted stage, each with nonempty value provenance, and
- an explicit external pointwise cutoff evaluator plus provenance.

Sparse sets such as stages `{1,3}` are rejected rather than being called a truncation prefix. Duplicate stages, mixed schedules, missing/extra value payloads, bad `q`, non-finite values, and cutoff values outside `[0,1]` fail closed.

For `a_j q >= 1`, the existing admitted paper cutoff-support theorem is sufficient to set that stage contribution to exact zero without invoking the external evaluator. For `a_j q < 1`, the evaluator is used numerically only. The implementation therefore records `paper_fixed_cutoff_pointwise_evaluator_verified=false` even when a caller supplies a smooth cutoff implementation.

The returned object exposes the stage-by-stage cutoff arguments, weights and contributions together with the accumulated finite-prefix `(A,B,p)` correction. It always preserves:

- `actual_correction_field_values_verified=false`,
- `paper_fixed_cutoff_pointwise_evaluator_verified=false`,
- `infinite_schedule_constructed_here=false`,
- `eq_9_21_infinite_sum_constructed=false`,
- `proposition_9_9_verified=false`, and
- `paper_exact_velocity_available=false`.

Thus caller-supplied values or a generic cutoff can be used to test plumbing but cannot be promoted to paper-exact data.

## Independent regression boundary

`tests/test_section9_finite_prefix.py` uses manufactured values explicitly labeled **not paper data** and an independent test-only linear cutoff that is explicitly **not the manuscript cutoff**. At `q=1/50` with admitted scales `(8,16,40,80)`, it computes the cutoff arguments independently as `(4/25, 8/25, 4/5, 8/5)` and checks the manual weights `(1,1,2/5,0)` and resulting sums. The final stage's test evaluator raises if called at or beyond support, so success independently confirms the `a_j q>=1` theorem short-circuit.

Additional regressions reject sparse and duplicate prefixes, inexact value coverage, mixed schedule roots, boundary/out-of-domain `q`, invalid cutoff output, malformed vectors, non-finite scalars and empty value provenance.

These tests validate finite-prefix assembly and fail-closed boundaries only. They do not prove that the manufactured values are manuscript corrections and do not establish any convergence statement.

## Remaining blocker

A genuine Eq. (9.21) construction still requires actual Section 7/8 correction fields tied to the admitted theorem witnesses, the manuscript's fixed cutoff evaluator, the hierarchy-derived infinite scale schedule on a common domain, arbitrary-order Eq. (9.20) flat residual bounds, and a proof that the locally finite infinite field sums are the Proposition 9.9 limit. Until those inputs are connected, this module must remain a finite formal-structure evaluator and `paper_exact_velocity_available` must remain false.
