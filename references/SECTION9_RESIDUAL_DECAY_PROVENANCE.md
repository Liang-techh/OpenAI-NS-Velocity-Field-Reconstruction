# Section 9 residual-decay ledger provenance

Status: **formal-structure only**.

This increment does not materialize the Section 9 correction sequence and does
not change `paper_exact_velocity_available=false`.

## Primary source

Paper: **Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026),
Section 9, especially Eqs. (9.17)–(9.21) and Proposition 9.9.  The repository's
canonical source URL is recorded in `references/SOURCES.md`.

The source snapshot used elsewhere in this repository records upstream
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd` as a
formalization cross-check.  No claim is made here that a Lean theorem directly
exports the executable arithmetic below.

## Paper identities encoded exactly

`section9_residual_decay.py` records, in rational arithmetic,

- `sigma_0 = 1/5` and `sigma_{j+1} = sigma_j + 1/10`, hence
  `sigma_j = 1/5 + j/10`;
- the Eq. (9.17) correction gain `g_j = h*j/10`;
- the admissible derivative loss
  `ell_m = 2*A + (m+1)*(1+3h/2)`, with `A=1/2+h`;
- the Eq. (9.18) leading residual exponent `h*sigma_j-K_m`; and
- exact ceiling formulas for the first stage where either exponent reaches a
  requested nonnegative power of `q`.

This is useful for auditing the eventual Section 9 implementation because the
paper's gains diverge linearly with `j`.  The code uses `Fraction` rather than
binary floating point for all schedule and threshold arithmetic.

## Independent regression checks

The tests substitute concrete rational data directly into the paper formulas,
without asking the production threshold helpers for an oracle.  In particular,
for `h=1/200` and derivative order `m=2`, they independently obtain
`ell_m=1613/400` and verify that `j=10065` is the first stage with
`g_j-ell_m >= 1`.  Separately, with `K_m=3/2`, they verify from
`h*(1/5+j/10)-K_m` that `j=6998` is exactly the first stage reaching residual
power `q^2`, while `j=6997` does not.

The returned audit record is explicitly fail-closed:
`status="formal-structure"`, `actual_residual_bound_verified=false`,
`flat_remainder_verified=false`, and `paper_exact_velocity_available=false`.

## Truth boundary / remaining blocker

The paper's estimates contain substantially more information than this ledger.
This increment does **not** establish:

- actual Section 7/8 corrections `Z_j` or `Delta u_j` satisfying Eq. (9.17);
- any constants `C_{j,m}`, logarithmic powers `P_{j,m}`, or derivative losses
  `K_m` from the genuine completed fields;
- the super-polynomial remainder bound
  `E_{j,m} <= C_{j,m,N} q^N` for every `N`;
- local finiteness and cutoff compatibility of the Eq. (9.21) potential sums;
- convergence of the actual velocity/pressure correction sequence; or
- Proposition 9.9 / Eq. (9.20) flatness for the repository's eventual residual.

In particular, `K_m` is deliberately caller-supplied as a theorem datum.  It is
never fitted from residual samples, and crossing a power threshold is not a
certificate that the physical residual has that decay.

The next high-value Section 9 step, once genuine Section 7/8 correction data are
available, is to attach measured/theorem-derived norm bounds to this exact
schedule and verify the hypotheses of (9.17)–(9.18) stage by stage before any
summation or convergence claim is enabled.
