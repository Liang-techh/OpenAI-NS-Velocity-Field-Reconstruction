# Primary amplitude ODE algebra provenance

Status: **formal-structure**. This file does not certify a paper-exact wave.

## Pinned sources

Official formalization: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

- `NavierStokes/MovingFrameODE.lean`
  - `coeff11`, `coeff12`, `coeff21`;
  - `rhsX`, `rhsY`;
  - `modal11`, `modal12`, `modal21`, `modal22`.
- `NavierStokes/GrowingMode.lean`
  - `modalOperator`.
- `NavierStokes/PrimaryODE.lean`
  - `FrameData.errorA/errorB/errorC`;
  - `FrameData.damping`, `FrameData.coefficient`;
  - `FrameData.forceX/forceY/forcing`;
  - `solution_hasDerivAt` and `ambientSolution_hasDerivAt`.

The official `PrimaryODE.lean` describes its ambient reconstruction as the exact reconstruction into equation (27), with projected physical forcing and harmonic-dependent scalar viscosity. The present Python increment implements only the pointwise finite-dimensional algebra appearing before the Volterra solution is constructed.

## Executable mapping

`src/openai_ns_reconstruction/primary_amplitude_ode.py` records the scalar projections needed by the official moving-frame reduction. For a supplied datum it computes

- `a = rho (g_K-rho_dot)/(1+rho^2)`,
- `b = (2 F N_theta-rho rotation)/(1+rho^2)`,
- `c = -(2 F N_theta+g_N)+rho rotation`,
- the reference-mode defects `B=b-lambda/h`, `C=c-lambda h`,
- the four official modal error entries,
- `d=j^2 viscosity`,
- the exact `GrowingMode.modalOperator` matrix,
- the projected forcing components and their moving-eigenbasis transform.

It also exposes the invertible coordinate change `x=p+q`, `y=h(p-q)` and its differentiated form using `h'=eigenRate*h`. This permits a regression to compare the modal ODE against the independently evaluated physical `rhsX/rhsY` equations rather than checking a matrix against itself.

## Independent regression boundary

`tests/test_primary_amplitude_ode.py` computes `rhsX/rhsY` directly from the physical moving-frame formulas and separately computes the modal RHS. It differentiates the basis change and checks that both paths agree. The test also checks the forcing transform separately and fails closed for a zero moving-basis eigenvector, non-finite inputs, and non-integer harmonic indices.

The numerical fixture is only an algebraic regression point. It is not a paper pulse, a fitted coefficient family, or evidence for any Section 7 estimate.

## What is still missing

This increment does **not**:

1. instantiate the datum from the paper-exact Proposition 5.5 background / certified phase-frame data;
2. construct or certify the true pulse forcing and order-zero target `T_{0,*}`;
3. solve the finite-interval Volterra equation used by `PrimaryODE.solution` or prove the smooth/weighted ODE estimates;
4. connect the resulting actual amplitude to the landed covariance coefficients and supported-curl realization;
5. prove support, zero-germ, cylindrical derivative identities, mean correction, or Section 9 residual improvement.

Accordingly Stage 3–6 remains **formal-structure** and `paper_exact_velocity_available` must remain `false` until those upstream and downstream obligations are discharged.
