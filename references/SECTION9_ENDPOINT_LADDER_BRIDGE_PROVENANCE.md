# Section 9 finite endpoint-ladder bridge provenance

Status: **formal-structure only**.

This increment does not construct the genuine Section 9 correction sequence, does not prove any uniform late-time residual derivative estimate, and does not change `paper_exact_velocity_available=false`.

## Primary sources

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026), especially Section 9.4 / Lemma 9.8 / Eq. (9.18), Proposition 9.9, and the Section 10 endpoint smooth-extension step. The canonical paper URL is recorded in `references/SOURCES.md`.

The repository also pins `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`. No claim is made here that the new finite ladder interface itself is a theorem already machine-checked in Lean.

## Landed prerequisites

PR #100 added `section9_endpoint_bridge.py`. For one endpoint jet `D^n R`, it requires a separately justified physical-time estimate

`sup_{x in K} ||partial_t D^n R_j(t,x)|| <= C (1-t)^(-alpha)`

on one late-time compact, with Section 9 derivative order exactly `n+1`, `0 <= alpha < 1`, official Section 10 plateau start `t >= 3/4`, and endpoint `T=1`. A pointwise Eq. (9.18) certificate cannot enter that bridge.

PR #104 added `section9_flat_remainder.py`. It records a coherent **finite** family of theorem inputs of the form

`E_{j,m}(q) <= C_{j,m,N} q^N`

on one small-`q` domain. Its provenance explicitly forbids extrapolating a finite prefix to the all-orders flatness required by Proposition 9.9.

## What this increment adds

`section9_endpoint_ladder_bridge.py` batches the PR #100 bridge across exactly the contiguous endpoint derivative degrees `0,...,N`.

Every row must be a `Section9UniformEndpointDerivativeWitness`. All rows must use the same Section 9 stage and the same compact spatial-window identifier. Each row is independently admitted through `admit_section9_uniform_endpoint_witness`; the resulting majorants are then packaged by the already-landed `Section10EndpointMajorantLadder`. Per-degree validity starts may differ, and the existing Section 10 ladder correctly takes their maximum as the common late-time start.

No coefficient, singularity exponent, validity start, or provenance string is fitted from residual samples by this module.

## Fail-closed q-flat / physical-time boundary

A `Section9FlatRemainderLadderCertificate` is rejected explicitly. Even a valid finite collection of `q`-flat estimates is not the theorem needed by the Section 10 Cauchy argument.

In particular, the repository still lacks the genuine common-domain/support result that would justify converting the Section 9 variables and residual estimates into a uniform physical-time estimate on the fixed Section 10 compact. The identity relating similarity variables (including the factor depending on `eta`) cannot be discarded to infer a global `q -> (1-t)` power law. PR #100 already records this boundary for one degree; this increment preserves it for a finite ladder.

## Regression boundary

`tests/test_section9_endpoint_ladder_bridge.py` checks that:

1. an out-of-order family is normalized only when it covers exactly endpoint degrees `0..N`;
2. the required Section 9 derivative orders remain `1..N+1` after single-degree admission;
3. mixed stages or spatial windows, missing powers, duplicates, and extra powers fail closed;
4. distinct admissible late-time starts are retained and the Section 10 common start is their maximum; and
5. a passing PR #104 finite `q`-flat ladder is still rejected as insufficient physical-time evidence.

These are interface-integrity tests, not evidence for the missing Section 9 theorem inputs.

## Remaining blocker

The repository still needs the actual Section 7/8 correction data and Section 9 residual sequence, the common shrinking-support/domain theorem, theorem-derived uniform late-time derivative bounds for every spacetime derivative on the relevant compact, locally uniform endpoint limits to all orders, genuine Borel gluing through `t=1`, smooth compact forcing, and the final bounded-energy / blow-up closure.

Accordingly this finite ladder bridge remains `formal-structure`, and `paper_exact_velocity_available` must remain `false`.
