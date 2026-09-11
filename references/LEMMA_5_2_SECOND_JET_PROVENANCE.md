# Lemma 5.2 repaired second-jet provenance

Status: **formal-structure only**.

This increment does not materialize the profile-dependent Section-5 hierarchy and
does not change `paper_exact_velocity_available=false`.

## Pinned source and landed dependencies

Official Lean source commit:
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Paper locations: Lemma 5.2, Eqs. (5.14)-(5.16), with downstream recursive use
through Eqs. (5.2)-(5.7).

The repository already had:

- `background_moment_repair.py` for the two compact `U_n` bumps, three compact
  `E_n` bumps, and the pointwise Eq. (5.16) solve;
- `background_moment_repair_jets.py` for arbitrary finite ordinary eta
  derivatives of the five repair coefficients from supplied unrepaired moment
  jets and `p(eta)=e_* f(eta)` jets;
- `background_moment_repair_profile.py` for actual repaired `U_n(X,eta)` and
  `E_n(X,eta)` functions; and
- `background_lower_history_source.py`, whose strict-lower PositiveAxis source
  consumes `ProfileSecondJet` data at orders below the current recursion order.

## What this increment adds

`background_moment_repair_second_jets.py` closes the derivative-interface gap
between the function-level compact repair and the strict-lower-history source.
For a normalized compact radial bump

`b(R) = exp(4 - 1/(s(1-s))) / N`, `s=(R-left)/(right-left)`,

it analytically evaluates `b`, `d_X b`, and `d_X^2 b` after the paper variable
change `R=sqrt(2X)`.  On the open support,

`d_X b = b_R / R`,

`d_X^2 b = b_RR/R^2 - b_R/R^3`.

Because every Lemma-5.2 bump has support strictly away from `R=0`, these chain
rules are nonsingular wherever the correction is active.  Outside and at the
flat support boundary the correction jet is returned as exact zero.

Combining these radial derivatives with the landed ordinary eta derivatives of
`alpha_j(eta)` and `beta_j(eta)` yields second `(X,eta)` jets of the repaired
`U_n` and `E_n`, including the mixed `d_X d_eta` and `d_eta^2` entries required
by downstream second-jet interfaces.

## Independent regression checks

The tests do not differentiate the production formulas again.  They evaluate
the repaired `U_n`/`E_n` functions through the existing function-level adapter
and independently recover all first, pure-second, and mixed derivatives with
centered finite-difference stencils.  They also verify that:

1. off the five compact supports the returned second jet is exactly the
   caller-supplied unrepaired jet;
2. no unavailable eta-repair derivatives are demanded where the compact
   correction is identically flat; and
3. malformed base-jet providers fail closed.

Finite differences are test oracles only; production derivatives are analytic.

## Truth boundary / remaining blocker

The new adapter still accepts the unrepaired `U_n`/`E_n` second jets and the
moment/patch eta jets from its caller.  It therefore does **not** prove that
those inputs arise from the converged Eq. (5.7) coefficient or from the
manuscript's analytic hierarchy.

It also deliberately does **not**:

- identify the repaired `E_n` jet with the PositiveAxis `phi_n` variable without
  the corresponding paper/Lean semantic bridge;
- construct the second jet of the regular flux `beta_n=V_n/X`, which requires
  differentiating the landed Eq. (5.2) reconstruction with sufficient higher
  eta data;
- provide hierarchy-derived five moment functions or a uniform nonvanishing
  certificate for `e_* f`;
- prove a common analytic strip, converged positive-order Picard series,
  infinite recursive cutoff schedule/local finiteness, or Proposition 5.3
  all-jets-flat residual decay.

The next high-value recursive step is to connect genuine solved/repaired
coefficient data to the PositiveAxis history semantics, in particular the
regular-flux second jet, without introducing sampled or caller-chosen
surrogates.  Until that chain is complete, Stage 2 remains `formal-structure`.
