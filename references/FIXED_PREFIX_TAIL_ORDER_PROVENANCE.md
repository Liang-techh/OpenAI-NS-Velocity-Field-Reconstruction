# Fixed-prefix tail-order certificate provenance

Status: **formal-structure**.

This increment maps the finite-prefix tail exponents already proved in the pinned official Lean repository at commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd` into a small executable arithmetic gate.  It does not construct the missing positive-order coefficients and does not claim Proposition 5.3 residual flatness.

## Official theorem map

Source: `NavierStokes/SlowBorelBase.lean`.

- `uncut_fixed_prefix_bound`: for an admissible scale sequence and actual smooth coefficient family, the difference between `slowSum` and the order-`J` uncut prefix satisfies, for every `m <= J+3`,

  `||D^m(slowSum - uncutPrefix_J)|| <= 2^(-J) q^(h(J+1)-m)`

  on a sufficiently small positive-`q` neighborhood.

- `powered_fixed_prefix_bound`: after the physical pullback and leading power `b`, a common derivative budget `M <= J+3` has exponent

  `h(J+1) + b - 2M`.

- `exists_ordinary_uncut_tail` and `exists_powered_physical_tail_finite`: by increasing `J` (with the theorem choosing `J >= M`) one can meet any prescribed finite decay exponent.

The executable module `src/openai_ns_reconstruction/background_tail_order.py` keeps these inequalities at the exact arithmetic level for the supplied binary64 parameters by using `Fraction.from_float`.  `minimal_prefix_for_chart_target` and `minimal_prefix_for_physical_target` solve the theorem-side exponent inequalities without a floating tolerance.  `certify_schedule_physical_tail` additionally fails closed unless the already-constructed finite `SlowBorelCutoffSchedule` actually reaches the selected prefix order.

## What this does certify

Given the *independently established* hypotheses of the Lean theorems, the module records the exact dyadic prefactor `2^-J`, the chart exponent `h(J+1)-m`, the powered physical exponent `h(J+1)+b-2M`, and the minimal finite prefix order needed for a requested exponent target.  Tests independently hand-check the rational inequalities and verify that a too-short finite schedule is rejected.

## What remains open

This increment does **not** certify any of the following:

- the actual profile-dependent `A0/A1/f_n` inputs or convergence of the Eq. (5.7) Picard solve;
- the materialized recursively repaired coefficient hierarchy;
- the true compactness bounds `C[j,m]` required to make the full scale sequence admissible;
- the infinite recursive schedule or theorem-level local finiteness;
- the finite coefficient identities that cancel the Navier--Stokes residual order by order;
- Proposition 5.3 / `BaseResidual` all-jets-flat residual decay.

Accordingly Stage 2 remains `formal-structure` and `paper_exact_velocity_available` must remain false.  The certificate is a reusable truncation-order bridge that becomes applicable to the actual residual only after the missing coefficient identities and admissibility hypotheses have been materialized.