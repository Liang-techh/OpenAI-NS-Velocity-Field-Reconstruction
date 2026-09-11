# Primary amplitude finite-interval solution provenance

Status: **formal-structure**. This file does not certify a paper-exact pulse or velocity field.

## Pinned source

Official formalization: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The relevant source is `NavierStokes/PrimaryODE.lean`:

- `extendedFamily` is the differentiable extension of the finite-interval Volterra solution;
- `extendedFamily_initial` fixes the initial value;
- `extendedFamily_hasDerivAt` identifies the derivative with `A(t) z(t) + f(t)` on the closed interval;
- `solution` specializes that construction to `FrameData.coefficient j` and `FrameData.forcing f`;
- `solution_initial` and `solution_hasDerivAt` record the corresponding IVP statements;
- `ambientSolution_hasDerivAt` reconstructs the modal solution into equation (27) once the official kinematic hypotheses hold.

The landed `primary_amplitude_ode.py` already transcribes the pointwise coefficient matrix, forcing transform, moving-basis reconstruction, and an independent physical/modal algebra check. The present increment addresses the next finite-interval interface only.

## Executable mapping

`src/openai_ns_reconstruction/primary_amplitude_solution.py` adds:

- `PrimaryAmplitudeInterval`, a fail-closed closed interval plus integer harmonic;
- `primary_modal_rhs_at`, which evaluates the landed pinned modal algebra against time-dependent `PrimaryModalDatum` and projected forcing callbacks;
- `solve_primary_amplitude_ivp`, a dense SciPy DOP853 numerical path for the resulting finite-interval IVP;
- `volterra_diagnostic`, which treats a supplied trajectory as a black box and independently evaluates the integral equation
  `z(t) - z(a) - integral_a^t (A(s) z(s) + g(s)) ds` with `quad_vec`.

The numerical solver and the quadrature diagnostic are intentionally separate. A trajectory is not accepted merely because it came from the solver: the Volterra checker recomputes the time integral from the callback path. It also requires the initial value exactly and fails closed on out-of-interval queries, malformed forcing, non-finite values, invalid tolerances, and non-integer harmonic labels.

## Independent regression

`tests/test_primary_amplitude_solution.py` uses a constant datum chosen so the official reference-mode subtraction makes the modal operator exactly `-d I`, with `d = j^2 viscosity`. For constant projected forcing the exact trajectory is then computed independently as

`z(t) = exp(-d tau) z0 + (1-exp(-d tau)) g/d`.

The test compares the numerical finite-interval path against that closed-form oracle at several times. It separately applies the Volterra diagnostic to the analytic trajectory, then perturbs one component by a quadratic term and checks that the integral defect becomes nonzero. Thus the regression is not a solver-output-versus-itself equality and does not use `f=R` followed by `R-f=0`.

The constant datum is only an algebraic regression specialization. It is **not** a hand-tuned surrogate for the manuscript pulse.

## Boundary / remaining obligations

This increment does **not** prove or claim that arbitrary callbacks satisfy the hypotheses of `PrimaryODE.solution`. In particular it does not establish:

1. that `datum_at` comes from the paper-exact Proposition 5.5 background, slow labels, phase, and moving tangent frame;
2. continuity/smoothness of the coefficient or forcing paths, the `FrameData.Kinematics` derivative identities, or Eqs. (7.9)-(7.11) uniform bounds;
3. the true pulse forcing, order-zero covariance target, or any paper-selected initial datum;
4. rigorous error bounds for SciPy `solve_ivp` or `quad_vec` (the reported quadrature error is a numerical estimate, not interval arithmetic or a Lean proof object);
5. supported-curl smoothness/support/zero-germ hypotheses, compact mean correction, or Section 9 residual-improvement convergence.

Accordingly this module is reusable **formal-structure** plumbing only. Stage 3-6 remains formal-structure and `paper_exact_velocity_available` must remain `false` until the genuine upstream data and all downstream certificates are materialized.
