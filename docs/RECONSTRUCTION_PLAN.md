# Reconstruction plan and truth-status ledger

The target is not merely to reproduce the blow-up exponent.  The target is an executable counterpart of the velocity field used in OpenAI's proof, with every implemented layer tagged by provenance and tested independently.


## 0.2.0 implementation update — partial, not completed

The historical stage definitions below remain the project goal. The machine-readable
current ledger is returned by `python -m openai_ns_reconstruction status`.

- Stage 0 now has a scale-relative root solver, direct tau APIs and analytic chain rules.
- Stage 1 has improved quadrature, pressure/axis support and pointwise Appendix A
  moment primitives. **The actual leading-profile constructor remains missing.**
- Stage 2 now has an analytic finite background velocity including the derivatives
  of the q-dependent cutoff; the coefficient recursion remains missing.
- Stage 7 now supports an analytic localization product rule and a C-infinity
  experimental spatial cutoff; no smooth-through-t=1 force is instantiated.
- Stage 8 includes independent manufactured solutions, step refinement, finite-cylinder
  energy, source/status reports, a completion gate and labelled toy-data exports.
- Stages 3–6 are still pending. No profile or whole construction is promoted to
  `paper-exact` by these engineering changes. No Lean compilation was performed.

The original phrase `paper-exact` below describes correspondence of the symbolic
formula to a paper equation, **not exact floating-point arithmetic or a completed
paper-instance certificate**. See `NUMERICS.md` and `NEXT_TASKS.md` for the precise
implementation boundary and next constructive dependencies.

## Status convention

- **paper-exact** — formula/algorithm is a direct implementation of the published construction.
- **formal-structure** — the published algebraic structure is encoded, but one or more recursively chosen profiles/cutoffs/frequencies are still inputs.
- **diagnostic-only** — a numerical device used for checking or visualization; it is not part of the proof construction.
- **toy** — illustrative data that must never be presented as OpenAI's velocity field.

## Stage 0 — similarity geometry — paper-exact

Implemented in `coordinates.py`.

From Eq. (4.1):

```text
A = 1/2 + h,  D = 1/2 - h,  tau = 1-t,
z = q^D eta,  tau = q(1-eta^2),  X = r^2/(2q).
```

The code solves the equivalent monotone scalar equation

```text
q - z^2 q^(2h) = tau
```

and checks both defining identities numerically.

## Stage 1 — leading axisymmetric velocity — kinematics paper-exact; profile data pending

Implemented in `profiles.py` and `velocity.py`.

From Eqs. (4.3), (4.6), (4.7):

```text
u_theta^(0) = q^(-A) E
u_z^(0)     = q^(-A) U
r u_r^(0)   = V0

A_X(U) = X^-1 integral_0^X U(x,eta) dx
V0 = X/L [2 eta U - 2D eta A_X(U) - (1-eta^2) d_eta A_X(U)]
L = 1 - 2h eta^2.
```

### Remaining work for a literal OpenAI leading profile

Implement the constructive content behind Theorem 4.6 rather than substituting an arbitrary `E,U` pair.  This requires translating the outer heat profile, inner analytic profile, matching moments, shear modification, cone constraints, and the Appendix A/B/C constructions into executable choices.  All numerical constants must be emitted to a provenance manifest.

## Stage 2 — all-order corrected background — formal-structure implemented

`background.py` encodes the Section 5 expansion with `lambda_n = 2nh` and the Stokes-streamfunction/vector-potential representation.  The exact recursive solver for `(phi_n,U_n,Pi_n,V_n)` from Eqs. (5.2)–(5.6), compact moment correction, and the paper's smooth cutoff sequence remain to be implemented.

Acceptance criterion: for every requested finite derivative order, the executable residual of the cutoff-summed background must show the predicted super-algebraic decay in `q`, and each coefficient must be reproducible from saved construction parameters.

## Stage 3 — dyadic charts, phases, and transported waves — pending

Translate Section 6 and the phase dynamics used by Section 7:

- dyadic scale `Q=2^-ell`, `epsilon=Q^h`, slow scale `S_* = ell^2`;
- chart coordinates and physical/normalized conversion;
- slow labels and separated auxiliary-torus supports;
- transported wavevectors/polarizations and rounded frequencies.

No arbitrary frequency list is acceptable for `paper-exact` status.

## Stage 4 — oscillatory realization of residual stress — pending

Translate Section 7, including:

- admissible stress cone;
- positive covariance representation;
- amplitude solve around the primary wave field;
- exact divergence-free realization by vector potentials/curls;
- curl remainder terms;
- physical evaluation map and support separation.

Acceptance criterion: averaged quadratic momentum flux reproduces the requested stress to the paper's stated order and all added wave velocities are divergence-free by construction.

## Stage 5 — compact mean corrections — pending

Translate Section 8 radial inverses, moment-preserving compact corrections, and the five-equation defect solve.

Acceptance criterion: zero auxiliary-average residual terms are cancelled while the support and moment conditions required by the later iteration remain exactly satisfied.

## Stage 6 — residual-improvement iteration — pending

Translate Section 9's full correction cycle and summation:

```text
Delta u_j = curl(A_j) + B_j
A = A0 + sum_j chi(a_j q) A_j
B e_theta = B0 e_theta + sum_j chi(a_j q) B_j
u_loc = curl(A) + B e_theta.                    (9.21)
```

The cutoff scales `a_j` must be chosen recursively from the paper's estimates, not manually tuned to plots.

Acceptance criterion: the local momentum residual and every tested derivative vanish to increasing order at the singularity, with convergence stable under increased truncation depth.

## Stage 7 — final compact field and force — composition implemented; paper data pending

`local_field.py` encodes the final structural operation

```text
u = curl(c A) + c B e_theta.                     (10.4)
```

The paper-specific spatial/time cutoff and the completed local fields must still be supplied.  Force reconstruction must use

```text
f = u_t + (u.grad)u - Delta u + grad p
```

at viscosity one, together with the smooth extension through `t=1` described in Section 10.

## Stage 8 — independent verification and visualization — started

`verify.py` provides numerical divergence and residual diagnostics.  Implemented numerical additions include analytic coordinate/background derivatives,
independent manufactured-solution tests, Taylor–Green refinement, finite-cylinder
energy, direct-tau exponent probes and labelled NPZ vortex slices. Still missing are
a general symbolic/automatic-differentiation backend, particle trajectories and a
reviewed cross-check table against the official Lean theorem modules.

## Non-negotiable rule

A toy profile may demonstrate the coordinate geometry, but it must never be called the OpenAI counterexample.  The repository reaches its stated goal only when all stages above are instantiated by the choices in the paper and can regenerate the final field from a clean checkout.
