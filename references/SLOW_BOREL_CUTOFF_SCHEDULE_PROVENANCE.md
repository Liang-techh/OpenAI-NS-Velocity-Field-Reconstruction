# Section 5 recursive cutoff-scale schedule provenance

## Source

Primary paper: *Finite Time Blowup for Navier-Stokes* (OpenAI, September 2026), Section 5 all-order background construction and the later diagonal/Borel cutoff argument used to sum the coefficient hierarchy.

Official formalization pinned by this repository: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant Lean modules:

- `NavierStokes/SlowBorelBase.lean`
- `NavierStokes/DiagonalScale.lean`
- downstream use in `NavierStokes/ConstructedSlowBase.lean`

`SlowBorelBase.AdmissibleScales` requires a positive integer schedule `a_j`, at least doubling at every step, and for each positive order `j` and derivative order `m <= j+2` both ordinary and blown-jet estimates with dyadic factor `2^(-j)`.

`SlowBorelBase.exists_admissibleScales` obtains normalized-template derivative bounds `C[j,m]` from compactness, then invokes `DiagonalScale.exists_diagonal_scales` with logarithmic exponent `p=0` and gain `g(j)=2*h*j`. In this specialization, the scalar smallness row reduces to

`C[j,m] q^(h*j) <= 2^(-j)` for `0 < q <= 1/a_j`.

`DiagonalScale.doublingEnvelope` is the recursive numerical schedule mechanism: a local stage scale is first chosen, then the final sequence dominates it and at least doubles from one stage to the next.

The same pinned `SlowBorelBase.lean` also proves the support side needed to turn the formal infinite sum into a locally finite object. `slowStage_locally_zero` obtains an eventual zero tail in a neighborhood of every positive-`q` chart point, and `slowSum_finite_at_scale` states that at fixed `q>0` the slow sum is represented by a finite prefix independently of the coefficient values. This uses the fixed smooth-cutoff geometry together with scales tending to infinity.

## Executable mapping

`src/openai_ns_reconstruction/background_cutoff_schedule.py` provides a finite-prefix executable witness for the scale-selection specialization:

- `local_scale_from_template_bounds(h, j, C_j)` takes the full triangular set `C[j,m]`, `m=0,...,j+2`, and chooses an integer local scale from the closed-form `p=0` edge inequality;
- `build_slow_borel_cutoff_schedule` applies the pinned doubling-envelope recurrence `a_0=max(1,B)`, `a_j=max(b_j,2*a_{j-1})` for the requested prefix;
- `build_slow_borel_cutoff_schedule_from_provider` gives a typed boundary for a future analytic provider of the actual compactness bounds;
- `edge_log_margin` certifies the scale inequality in log space without underflow; and
- `majorant_at` evaluates only the resulting scalar bound on its permitted active interval.

`src/openai_ns_reconstruction/background_cutoff_support.py` adds the finite-prefix support/plateau certificate implied by the same scale sequence and the repository's fixed cutoff (`chi=1` on `s<=1/2`, `chi=0` on `s>=1`):

- `classify_cutoff_support(schedule,q)` uses the exact integer scales and `Fraction.from_float(q)` to classify the supplied binary64 `q` value without a tolerance;
- `plateau_through` records the exact uncut positive-order prefix with `a_j q <= 1/2`;
- because `a_{j+1} >= 2 a_j`, there can be at most one order with `1/2 < a_j q < 1`;
- `zero_from` records the first constructed order with `a_j q >= 1`, after which every later order in the constructed prefix has identically zero cutoff; and
- `stable_truncation_order_within_prefix` records the smallest truncation whose value is unchanged by appending any later stages already present in that finite schedule.

The certificate deliberately refuses to infer an infinite zero tail when the requested finite prefix ends before `a_j q >= 1`. This mirrors the theorem's support logic while keeping the executable claim no stronger than the data actually constructed.

`tests/test_background_cutoff_schedule.py` independently recomputes the `p=0` edge inequality, checks a hand-solvable threshold sequence, verifies the recursive doubling envelope, verifies the exact triangular `(j,m)` provider calls, and exercises fail-closed malformed/forged schedules.

`tests/test_background_cutoff_support.py` independently evaluates the actual `standard_cutoff` on a hand-solvable doubling schedule, checks the closed `1/2` plateau and `1` zero boundaries, verifies the unique transition-order consequence of doubling, checks finite-prefix truncation stability once a zero tail is reached, and verifies that a prefix ending in the transition collar remains explicitly uncertified beyond that prefix.

## Truth boundary

This remains **formal-structure / solver infrastructure**, not the completed paper-exact all-order cutoff sequence.

The actual constants `C[j,m]` must come from uniform compactness bounds for the materialized recursively repaired Section 5 coefficient fields. The scale constructor deliberately does not estimate them from a finite sample grid or accept a generic cutoff as evidence. The current executable object materializes only a requested finite prefix; the support certificate proves exact plateau/transition/zero facts only inside that prefix. It does **not** establish the paper's infinite schedule, global local-finiteness theorem, arbitrary-order residual flatness, or Proposition 5.3 by itself.

Once the true coefficient hierarchy is available, its analytic normalized-template bounds can feed the provider interface and the same recurrence can instantiate longer paper-admissible prefixes; an actual infinite schedule plus the theorem-side convergence argument is still required before claiming full local finiteness/residual flatness. Until then Stage 2 remains `formal-structure`, and `paper_exact_velocity_available` remains false.
