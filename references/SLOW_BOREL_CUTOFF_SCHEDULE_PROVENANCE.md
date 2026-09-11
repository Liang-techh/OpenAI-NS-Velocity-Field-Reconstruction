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

## Executable mapping

`src/openai_ns_reconstruction/background_cutoff_schedule.py` provides a finite-prefix executable witness for that exact specialization:

- `local_scale_from_template_bounds(h, j, C_j)` takes the full triangular set `C[j,m]`, `m=0,...,j+2`, and chooses an integer local scale from the closed-form `p=0` edge inequality;
- `build_slow_borel_cutoff_schedule` applies the pinned doubling-envelope recurrence `a_0=max(1,B)`, `a_j=max(b_j,2*a_{j-1})` for the requested prefix;
- `build_slow_borel_cutoff_schedule_from_provider` gives a typed boundary for a future analytic provider of the actual compactness bounds;
- `edge_log_margin` certifies the scale inequality in log space without underflow; and
- `majorant_at` evaluates only the resulting scalar bound on its permitted active interval.

`tests/test_background_cutoff_schedule.py` independently recomputes the `p=0` edge inequality, checks a hand-solvable threshold sequence, verifies the recursive doubling envelope, verifies the exact triangular `(j,m)` provider calls, and exercises fail-closed malformed/forged schedules.

## Truth boundary

This is **formal-structure / solver infrastructure**, not the completed paper-exact all-order cutoff sequence.

The actual constants `C[j,m]` must come from uniform compactness bounds for the materialized recursively repaired Section 5 coefficient fields. This module deliberately does not estimate them from a finite sample grid or accept a generic cutoff as evidence. The current executable object materializes only a requested finite prefix; it does not claim the infinite schedule, local finiteness of the full series, arbitrary-order residual flatness, or Proposition 5.3 by itself.

Once the true coefficient hierarchy is available, its analytic normalized-template bounds can feed the provider interface and the same recurrence can instantiate the paper-admissible scales. Until then Stage 2 remains `formal-structure`, and `paper_exact_velocity_available` remains false.
