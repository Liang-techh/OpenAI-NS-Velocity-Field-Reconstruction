# Section 9 pointwise stage-certificate provenance

Status: **formal-structure only**.

This increment does not construct the Section 9 correction sequence and does
not change `paper_exact_velocity_available=false`.

## Primary source

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026),
Section 9.4, Lemma 9.8, printed pages 113–114, especially Eqs. (9.17)–(9.19).
The canonical paper URL is recorded in `references/SOURCES.md`.

The paper states, for every Cartesian space-time derivative order `m`,

- Eq. (9.17):
  `|Z_j|_m + |Delta u_j|_m <= C_{j,m} q^(g_j-ell_m)
   (1+|log q|)^(P_{j,m})`, with `g_j=h j/10`;
- Eq. (9.18):
  `|R(u^[j],p^[j])|_m <= C_{j,m} q^(h sigma_j-K_m)
   (1+|log q|)^(P_{j,m}) + E_{j,m}`;
- the same equation separately requires
  `E_{j,m} <= C_{j,m,N} q^N` for every `N`;
- the proof gives the admissible common derivative loss
  `ell_m = 2A + (m+1)(1+3h/2)` used by the existing exact exponent ledger.

The source snapshot used by the repository is
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
This increment does not claim that a Lean theorem directly exports these
pointwise numerical envelopes.

## What this increment adds

`section9_stage_certificate.py` is a fail-closed arithmetic adapter from a
**separately certified** scalar norm bound to the exact Section 9 exponent
schedule already implemented in `section9_residual_decay.py`.

For Eq. (9.17), the caller must supply a certified upper bound for
`|Z_j|_m + |Delta u_j|_m`, together with `C_{j,m}`, `P_{j,m}`, `q`, `h`, and
provenance.  The adapter computes the exact rational exponent
`g_j-ell_m` and checks the paper right-hand side.

For Eq. (9.18), the caller separately supplies a certified residual upper bound
and a certified pointwise upper bound for `E_{j,m}`.  The adapter computes the
exact exponent `h sigma_j-K_m`, keeps the flat remainder as a separate term,
and checks the sum of the two right-hand-side contributions.

The power comparison is done in logarithmic scale.  This avoids binary64
underflow for very large positive exponents without changing any theorem
constant or exact rational exponent.

`CertifiedBoundDatum` accepts only explicitly classified evidence kinds:
`paper-derived`, `lean-derived`, `certified-numerical`, or
`rigorous-external`.  Ordinary sampled/fitted values are rejected by the API.
A nonempty provenance string is mandatory.

## Independent regression checks

`tests/test_section9_stage_certificate.py` computes the paper right-hand sides
independently in ordinary scale for two closed arithmetic cases:

1. `h=1/200`, `m=2`, `j=10065`, where direct substitution gives
   `g_j-ell_m=1`; the test places certified data just below and just above the
   independently calculated Eq. (9.17) envelope.
2. `h=1/200`, `K_m=3/2`, `j=6998`, where direct substitution gives
   `h sigma_j-K_m=2`; the test independently adds the supplied flat remainder
   and again checks both sides of the threshold.

A separate regression uses an exponent large enough that direct `q**exponent`
would underflow, demonstrating that the production comparison stays in log
scale.  Invalid `q`, constants, derivative losses, evidence kinds, and missing
provenance fail closed.

## Truth boundary / remaining blocker

Passing this certificate is **not** a proof of Lemma 9.8.  In particular, this
increment does not:

- construct genuine Section 7/8 corrections `Z_j` or `Delta u_j`;
- derive `C_{j,m}`, `P_{j,m}`, or `K_m` from those fields;
- derive a certified norm bound from the field itself;
- prove the Eq. (9.17) or Eq. (9.18) envelope uniformly for all sufficiently
  small `q` or uniformly over the common Section 9 domain;
- prove `E_{j,m} <= C_{j,m,N} q^N` for every `N` from a single pointwise
  remainder datum;
- establish the common-domain statement of Lemma 9.7, the locally finite
  shrinking-cutoff sum (9.21), convergence of the actual correction sequence,
  or Proposition 9.9 flatness.

The certificate therefore intentionally reports
`uniform_in_q_verified=false`, `flat_remainder_all_orders_verified=false`,
`actual_section9_sequence_verified=false`, and
`paper_exact_velocity_available=false` even when its arithmetic comparison
passes.

The next promotion requires genuine upstream wave/mean correction data to
produce certified bounds over the whole common domain, followed by an
all-orders flat-remainder witness rather than pointwise data.
