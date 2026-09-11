# Section 9 Eq. (9.21) local-sum admission provenance

Status: **formal-structure only**.

This increment does not construct the Section 9 correction fields, does not prove Proposition 9.9, and does not change `paper_exact_velocity_available=false`.

## Primary source

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026), Lemma 5.4 and its local-finiteness proof, plus Section 9 / Proposition 9.9 Step 2 and Eq. (9.21). The canonical paper URL is recorded in `references/SOURCES.md`.

The paper forms

`A = A0 + sum_{j>=1} chi(a_j q) A_j`,

`B e_theta = B0 e_theta + sum_{j>=1} chi(a_j q) B_j`,

and the analogous pressure sum. The common small-`q` domain is supplied by the Section 9 construction, while Lemma 5.4 uses the recursive scale condition `a_{j+1} >= 2 a_j` and the fixed cutoff support rule `chi(s)=0` for `s>=1`. On a compact strip `c <= q <= c' < q_big`, `q` is bounded away from zero, so sufficiently large scales force every later cutoff factor to vanish.

## What this increment adds

`section9_local_sum_admission.py` encodes exactly that arithmetic implication without importing or reimplementing the Section 5 cutoff-scale constructor. `Section9LocalSumHypothesis` is an external theorem/proof datum: it asserts one common domain `0<q<q_big` and an **infinite** scale sequence whose first scale is `a_1` and whose successive scales at least double. The datum requires rigorous provenance metadata and rejects sampled/fitted evidence labels.

For a requested compact strip `0<c<c'<q_big`, `certify_eq_9_21_local_finiteness` uses exact `Fraction` arithmetic and the theorem inequality

`a_j >= a_1 2^(j-1)`.

A positive correction stage can remain active only if

`a_1 2^(j-1) c < 1`.

The certificate therefore computes the exact worst-case finite stage bound, including the strict cutoff edge where `a_j c = 1` is already inactive. It never infers an infinite schedule from a finite scale prefix and never inspects correction samples.

The returned certificate deliberately keeps

- `source_theorems_machine_verified=false`,
- `actual_correction_fields_verified=false`,
- `eq_9_21_sum_constructed=false`,
- `proposition_9_9_verified=false`, and
- `paper_exact_velocity_available=false`.

## Regression boundary

`tests/test_section9_local_sum_admission.py` independently enumerates the worst-case exact doubling sequence and compares it with the production closed arithmetic bound. It checks a generic finite-activity case, the exact support edge `a_j q=1`, a strip on which no positive correction can be active, and fail-closed rejection of invalid common domains, non-strict compact strips, sampled evidence, missing provenance, and invalid first-scale data.

These regressions verify only the local-finiteness implication from the admitted theorem hypotheses. They do not verify the external common-domain or infinite-schedule theorem itself.

## Remaining blocker

The repository still needs the genuine Section 7/8 correction fields `A_j`, `B_j`, `p_j`, their common-domain and smooth zero-extension/support hypotheses, the actual infinite scale schedule arising from theorem-derived stage bounds, the arbitrary-order flat residual conclusion of Eq. (9.20), and the resulting locally finite smooth sums of Eq. (9.21). Only after those inputs are real can this admission arithmetic participate in a genuine Proposition 9.9 convergence record.
