# Section 9 admitted-stage physical support provenance

Status: **formal-structure only**.

This increment does not construct the Section 7/8 correction values, does not construct the Eq. (9.21) sum, does not prove endpoint convergence, and does not change `paper_exact_velocity_available=false`.

## Primary source

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026), Eq. (4.1), Lemma 5.4, and Section 9 / Proposition 9.9, especially Eq. (9.21). The canonical source is recorded in `references/SOURCES.md`.

The already-landed Section 9 admission layers establish two separate formal facts:

1. `section9_correction_extension_admission.py` admits one positive correction stage only after theorem/proof witnesses cover the actual stage scale `a_j`, the common `q_big`, fixed cutoff support, and smooth zero extension for all three Eq. (9.21) components `A_j/B_j/p_j`.
2. `section9_physical_window.py` uses Eq. (4.1) and the official Section 10 support `|z|<=1/4` to enclose every point of a closed late-time pre-endpoint physical slab in an exact rational q-strip `q_lower <= q <= q_upper`.

Eq. (9.21) multiplies each positive-stage correction by `chi(a_j q)`, while the paper cutoff support rule gives `chi(s)=0` for `s>=1`. Therefore an admitted stage is identically zero on the entire fixed Section 10 support throughout a certified slab whenever

`a_j * q_lower >= 1`.

## What this increment adds

`section9_stage_physical_support.py` composes the two landed certificates above without sampling correction values. It first verifies that the admitted stage and the physical slab belong to the same common `q_big` and the same schedule first scale. It recomputes the exact doubling lower bound `a_1 2^(j-1)` and refuses a stage certificate derived from a different schedule.

It then computes the exact cutoff-argument interval

`a_j q_lower <= a_j q <= a_j q_upper`.

If the lower endpoint is at least one, the returned certificate records that the complete admitted `A_j/B_j/p_j` stage is uniformly zero on that physical slab. Failure of this inequality is **not** interpreted as positive stage activity.

Because the bridge uses the witnessed actual `a_j`, it can sharpen the generic local-finiteness bound that knows only `a_j >= a_1 2^(j-1)`. The implementation also cross-checks that whenever the generic doubling bound excludes a stage, the actual-scale support test excludes it as well.

## Independent regression boundary

`tests/test_section9_stage_physical_support.py` contains an exact example with `a_1>=8`, an admitted actual scale `a_1=12`, and the physical slab `7/8 <= t <= 11/12`. The Eq. (4.1) slab bridge gives `q_lower=1/12`, hence the actual stage satisfies `12*q>=1` everywhere and is uniformly zero, while the generic lower bound `8*q_lower=2/3` is insufficient to exclude it. This demonstrates a real strengthening from the admitted actual schedule value rather than a caller-selected surrogate.

A second regression independently calls the existing physical root solver for Eq. (4.1) on multiple `h,t,z` points and confirms the certified support inequality. Additional tests reject mixed common domains, corrupted schedule provenance, endpoint slabs, and any attempt to infer activity when uniform zero cannot be proved.

These tests validate the bridge and exact support arithmetic only. They do not machine-prove the external theorem witnesses used to admit a stage and do not synthesize correction values.

## Relation to Section 10

This is a paper-shaped time-localization fact, not an arbitrary user time window: the stage cutoff is exactly the Eq. (9.21) factor `chi(a_j q)`, and the physical-time enclosure comes from Eq. (4.1) plus the fixed Section 10 support geometry. It provides a safe zero short-circuit on closed pre-endpoint slabs and a stronger local-finiteness record for future construction of the Eq. (9.21) sum.

The bridge deliberately remains fail-closed at `t=1`. At `z=0`, Eq. (4.1) gives `q=1-t`, so no positive uniform `q_lower` survives to the endpoint. Consequently this increment does not provide the arbitrary-order physical-time residual derivative majorants, locally uniform endpoint limits, Borel glue, smooth compact forcing, finite-energy bound, or blow-up closure required by Issue #4.

## Remaining blocker

A genuine Section 10 completion still requires actual Section 7/8 correction fields and theorem-connected stage witnesses, construction and convergence of the locally finite Eq. (9.21) fields, the arbitrary-order Eq. (9.20) residual flatness/common-domain theorem, uniform physical-time derivative majorants up to `t=1`, and the resulting smooth forcing / energy / blow-up certificates. Until those are present, `paper_exact_velocity_available` must remain false.
