# Section 5 finite residual truncation provenance

Status: **formal-structure**.

This increment makes the finite-order residual bookkeeping in the pinned official Lean formalization executable.  It is a cross-validation bridge between the already landed coefficient-recursion infrastructure and the fixed-prefix tail-order machinery; it does **not** supply the missing profile-dependent coefficient hierarchy and does not claim Proposition 5.3 residual flatness.

## Official theorem map

Pinned official Lean source: `openai/NavierStokesAndEuler` commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, file `NavierStokes/SlowExpansionResidual.lean`.

The executable module `src/openai_ns_reconstruction/background_truncation_residual.py` mirrors the following definitions/theorems:

- `slowOrder`: `lambda_n = 2 n h`;
- `convolution` / `finiteConvolution_eq`: retained Cauchy-product coefficients;
- `pairTail` / `pair_sum_split`: the interactions with `i,j <= N` but `i+j > N` that remain after truncating the coefficient equations at order `N`;
- `previous` / `shifted_sum`: the one-order shift created by the axial-viscosity term, including the explicit final term at `N+1`;
- `recurrence`: `R_n = L_n + sum_{i+j=n} K_ij - previous(A)_n`;
- `recurrence_truncation`: the exact identity

  `finite residual = sum_{n<=N} w_n R_n + pairTail_N - w_{N+1} A_N`;

- `recurrence_truncation_of_zero`: once the actual retained coefficient equations are independently proved to vanish, the finite residual is exactly the omitted pair tail minus the final shifted-viscosity term.

The official file explicitly states that finite products are expanded exactly and that terms above the retained slow order and the last axial-viscosity term are displayed as finite remainders.

## Executable increment

`recurrence_truncation_breakdown(...)` evaluates both sides of the theorem without dropping a nonzero retained recurrence.  It reports the individual recurrence values, their weighted contribution, the omitted pair tail and the final shifted term.  There is deliberately no tolerance-based switch that silently turns an approximately small recurrence into the theorem hypothesis `R_n = 0`.

`omitted_slow_tail_majorant(...)` adds only the elementary finite-order consequence for the paper weights

`w_n = q^(b + 2 n h)`,  `0 < q <= 1`:

`|pairTail_N - w_{N+1} A_N|`
`<= q^(b + 2 (N+1) h) * (sum_{i+j>N} |K_ij| + |A_N|)`.

Thus, **conditional on separately established exact recurrence cancellation**, the finite residual has at least the first omitted slow-order factor.  The exponent is stored as an exact `Fraction` for the supplied binary64 `h` and `b`; numerical evaluation fails closed on binary64 underflow/overflow rather than returning a misleading zero or infinity.

Tests use an algebraic fixture whose `L_n` entries are independently chosen to make `R_n=0`, compare the production split against a separate direct finite expansion, verify that a deliberately nonzero retained recurrence is not discarded, and check the first-omitted-order majorant over several `q` values.  These fixtures are test algebra only, not paper coefficient data.

## Truth boundary

This increment does **not** establish any of the following:

- the actual profile-derived `A0/A1/f_n` or the converged Lemma 5.1 coefficients;
- that Eqs. (5.3)-(5.6) vanish for a materialized paper coefficient hierarchy;
- the true eta-dependent Lemma 5.2 repaired hierarchy and its support/stress conclusions;
- hierarchy-derived normalized jet bounds `C[j,m]` or the full infinite admissible cutoff schedule;
- derivative-level/chart-level estimates for the Navier--Stokes residual;
- Proposition 5.3 / `BaseResidual` all-jets-flat decay;
- a paper-exact velocity or force.

Accordingly Stage 2 remains **formal-structure** and `paper_exact_velocity_available` must remain `false`.  The new module becomes a paper-relevant residual cross-check only when the missing actual coefficient identities are fed into it rather than caller-supplied fixtures.
