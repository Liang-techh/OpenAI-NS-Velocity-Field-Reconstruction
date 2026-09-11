# Section 9 finite flat-remainder ladder provenance

Status: **formal-structure only**.

This increment does not construct the Section 9 correction sequence, does not prove the all-orders flat remainder, and does not change `paper_exact_velocity_available=false`.

## Primary source

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026), Section 9.4, Lemma 9.8, printed pages 113–114, especially Eq. (9.18). The canonical paper URL is recorded in `references/SOURCES.md`.

The paper's residual estimate separates the leading power-law term from `E_{j,m}` and requires, for every nonnegative integer `N`, a bound of the form

`E_{j,m}(q) <= C_{j,m,N} q^N`

on one sufficiently small common `q` domain. The source snapshot used by the repository is `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

PR #99 deliberately kept only a pointwise `E_{j,m}` upper bound in `section9_stage_certificate.py`; its provenance explicitly records that a pointwise datum is insufficient to establish all-orders flatness. PR #100 likewise refuses to promote a pointwise Eq. (9.18) certificate into a Section 10 endpoint theorem.

## What this increment adds

`section9_flat_remainder.py` adds the next fail-closed interface: an externally justified `Section9FlatRemainderPowerWitness` for a single power `N` must identify one Section 9 stage `j`, one derivative order `m`, one full small-`q` domain `0 < q <= q_upper`, a certified constant `C_{j,m,N}`, and theorem/proof provenance.

The witness is semantic theorem input: it means the same `E_{j,m}` satisfies `E_{j,m}(q) <= C_{j,m,N} q^N` on the whole stated domain. The module does not infer that statement from samples, finite differences, fitted slopes, or one supplied `q`.

`certify_flat_remainder_power_ladder` accepts a finite family only when it covers **exactly** the contiguous powers `N=0,...,Nmax`, with no duplicates or gaps, and when every row uses the same stage, derivative order and `q_upper`. This prevents a later integration step from silently combining unrelated domains or different residual stages and calling the result an all-orders flatness proof.

The returned `Section9FlatRemainderLadderCertificate` records only a coherent finite prefix. It therefore keeps

- `source_theorems_machine_verified=false`,
- `flat_remainder_all_orders_verified=false`,
- `actual_section9_sequence_verified=false`, and
- `paper_exact_velocity_available=false`.

## Regression boundary

`tests/test_section9_flat_remainder.py` checks the fail-closed structure independently of any residual sample:

1. an out-of-order but complete family `N=0,1,2,3` is normalized without inventing or dropping powers and preserves every supplied certified constant exactly;
2. missing, duplicate, or extra powers are rejected;
3. mixing Section 9 stages, derivative orders, or small-`q` domains is rejected; and
4. sampled evidence, missing provenance, invalid `q_upper`, and boolean pseudo-integers are rejected through the existing rigorous evidence gate.

These tests verify interface integrity only; they are not evidence that the external power-law theorems exist for the genuine Section 9 remainder.

## Remaining blocker

The repository still needs genuine Section 7/8 corrections and the actual Section 9 residual sequence to produce theorem-derived `E_{j,m}` power-law witnesses on the common domain for arbitrary `N`. A finite ladder cannot be extrapolated to all powers. The common-domain statement of Lemma 9.7, locally finite shrinking-cutoff sum (9.21), convergence of the actual correction sequence, and Proposition 9.9 flatness remain unresolved.
