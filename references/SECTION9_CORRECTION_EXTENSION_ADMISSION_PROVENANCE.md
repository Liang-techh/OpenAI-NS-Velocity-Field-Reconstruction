# Section 9 correction-extension admission provenance

Status: **formal-structure only**.

This increment does not construct the Section 7/8 correction fields, does not construct the Eq. (9.21) sum, does not prove Proposition 9.9, and does not change `paper_exact_velocity_available=false`.

## Primary source

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026), Lemma 5.4 and Section 9 / Proposition 9.9, especially Eqs. (9.20)–(9.21). The canonical paper URL is recorded in `references/SOURCES.md`.

Eq. (9.21) assembles the positive-stage vector-potential, swirl and pressure corrections with the fixed factors `chi(a_j q)`. The already-landed `section9_local_sum_admission.py` certifies only the exact local-finiteness arithmetic implied by a common small-`q` domain and an infinite schedule satisfying `a_{j+1} >= 2 a_j`. That arithmetic does not itself certify that a supplied stage correction has the smooth support/zero-extension properties required before it can participate in the locally finite field sum.

## What this increment adds

`section9_correction_extension_admission.py` adds a fail-closed theorem-input gate for one positive correction stage. A stage is admitted only after separate witnesses cover exactly the three Eq. (9.21) components `A_j`, `B_j`, and `p_j`. The witnesses must agree on:

- the positive stage `j`,
- the common domain `0 < q < q_big`,
- the stage cutoff scale `a_j`,
- rigorous provenance,
- membership of that scale in the admitted schedule,
- the paper cutoff-support statement, and
- smooth zero extension of the corresponding correction component.

The verifier then checks the machine-checkable consequence of the infinite doubling hypothesis exactly:

`a_j >= a_1 2^(j-1)`.

The schedule-membership assertion is kept separate from this inequality. Therefore an arbitrary caller-selected integer that merely satisfies the lower bound cannot be promoted to a paper stage unless an independent theorem/proof witness also certifies that it is the actual `a_j`.

The returned certificate deliberately keeps:

- `source_theorems_machine_verified=false`,
- `actual_correction_field_values_verified=false`,
- `infinite_schedule_constructed_here=false`,
- `eq_9_21_sum_constructed=false`,
- `proposition_9_9_verified=false`, and
- `paper_exact_velocity_available=false`.

## Independent regression boundary

`tests/test_section9_correction_extension_admission.py` checks the stage-3 schedule consequence independently as `8 * 2^(3-1) = 32` and verifies that an admitted scale `40` passes while scale `31` fails. It also verifies fail-closed rejection of missing or duplicate `A/B/p` coverage, inconsistent stage/domain/scale data, sampled evidence, empty provenance, and any unverified common-domain, schedule-membership, cutoff-support, or smooth-zero-extension fact.

These tests verify the admission logic and exact schedule arithmetic only. They do not prove the external theorem witnesses and do not inspect or synthesize correction values.

## Relation to upstream Lean artifacts

The pinned OpenAI Lean repository contains executable smooth-cutoff and wave-regularity infrastructure, including `ActualSignedStageControls.lean` and `ActualWaveRegularityData.lean`. Those files corroborate the proof architecture around smooth supported wave controls, but this increment does **not** treat them as a machine proof of the Section 9 `A_j/B_j/p_j` correction-extension hypotheses. Those hypotheses remain explicit external theorem inputs until their actual provenance is connected.

## Remaining blocker

The repository still needs genuine Section 7/8 correction fields and theorem-derived stage scales, a machine-connected common-domain/smooth-extension proof for those actual fields, the arbitrary-order Eq. (9.20) flat residual, and construction of the locally finite Eq. (9.21) sums. Only then can the existing local-finiteness, physical-window, and endpoint bridges participate in a genuine Proposition 9.9 convergence record.
