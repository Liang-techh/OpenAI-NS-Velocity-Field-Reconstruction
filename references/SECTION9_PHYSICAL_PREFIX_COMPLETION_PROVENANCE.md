# Section 9 physical-prefix completion provenance

Status: **formal-structure only**.

This increment does not construct the actual Section 7/8 correction fields, does not certify the executable pointwise cutoff as the manuscript's fixed cutoff, does not reach `t=1`, and does not change `paper_exact_velocity_available=false`.

## Primary source

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026), Lemma 5.4 and Section 9 / Proposition 9.9, especially Eqs. (9.20)–(9.21), together with the Eq. (4.1) similarity relation used by the landed Section 10 physical-slab bridge. The canonical source is recorded in `references/SOURCES.md`.

## Reused landed facts

This increment deliberately composes existing modules rather than reimplementing them:

1. `section9_local_sum_admission.py` checks the infinite schedule hypothesis `a_{j+1} >= 2 a_j` and derives the exact maximum stage that can remain active on a compact positive q-strip.
2. `section9_physical_window.py` transfers the official Section 10 support cylinder on a closed pre-endpoint slab `t0 <= t <= t1 < 1` into the analytic enclosure `1-t1 <= q <= 1-t0+1/16` and then reuses the local-sum admission.
3. `section9_correction_extension_admission.py` admits a concrete positive stage only after complete A/B/p common-domain, schedule-membership, cutoff-support, and smooth-zero-extension witnesses.
4. `section9_finite_prefix.py` evaluates a contiguous admitted finite prefix at one q while explicitly refusing to promote caller-supplied field values or a generic cutoff callable to paper-exact data.

## What this increment adds

`section9_physical_prefix_completion.py` closes one structural gap between those layers.

Given one physical-slab local-finiteness certificate and one already-evaluated contiguous prefix, it checks:

- the slab certificate came from the same `q_big`, `first_scale`, evidence kind, and evidence provenance as the supplied infinite-schedule hypothesis;
- the evaluation q lies inside the slab's certified analytic q enclosure;
- the supplied stage admission certificates cover exactly the evaluated prefix;
- every evaluated scale matches the corresponding admitted schedule-member scale; and
- the prefix reaches at least the slab certificate's `max_potentially_active_stage`.

If these checks pass, then every omitted positive stage begins no earlier than the already-certified `first_guaranteed_inactive_stage`. Therefore the paper cutoff support rule forces every omitted stage to zero on the **whole closed pre-endpoint physical slab**. This is a support-completeness statement for the stage index set; it is stronger than a pointwise truncation check but weaker than constructing the manuscript's Eq. (9.21) field.

The returned certificate keeps all critical truth gates false:

- `endpoint_covered=false`,
- `actual_correction_field_values_verified=false`,
- `paper_fixed_cutoff_pointwise_evaluator_verified=false`,
- `eq_9_21_infinite_sum_constructed=false`,
- `endpoint_uniform_residual_majorants_verified=false`,
- `proposition_9_9_verified=false`, and
- `paper_exact_velocity_available=false`.

## Independent regression boundary

`tests/test_section9_physical_prefix_completion.py` uses manufactured correction payloads explicitly labeled non-paper data.

For the slab `3/4 <= t <= 15/16`, the existing Eq. (4.1) bridge gives `q_lower=1/16`. With the independently supplied infinite-schedule root `a_1=8`, the local-finiteness arithmetic permits only stage 1 to remain active. The regression evaluates stage 1 at `q=1/10`, then independently checks for several omitted indices that

`8 * 2^(j-1) * (1/16) >= 1` for every `j >= 2`,

so the omitted tail is support-zero throughout the slab.

A second regression moves the slab closer to the endpoint, `7/8 <= t <= 127/128`, where `q_lower=1/128` and four stages may be active. A one-stage prefix is then rejected rather than mislabeled complete. Further tests reject an evaluation q outside the physical enclosure, a mismatched admitted scale, and a slab derived from different schedule provenance.

These tests validate arithmetic and provenance composition only. They do not prove that the manufactured payload is a manuscript correction field.

## Remaining blocker

A genuine Section 10 completion still requires machine-connected actual Section 7/8 correction values, the manuscript fixed cutoff evaluator (or an equivalent proved realization), the actual locally finite Eq. (9.21) fields, arbitrary-order residual derivative bounds uniform up to `t=1`, smooth endpoint extension, smooth compact forcing, bounded kinetic energy, and the final blow-up closure. Until those inputs exist, this bridge must remain formal-structure only.
