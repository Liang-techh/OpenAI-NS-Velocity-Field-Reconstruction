# OpenAI NS Velocity Field Reconstruction

Executable reconstruction of the finite-time Navier–Stokes blow-up velocity field described in OpenAI's September 2026 paper **Finite Time Blowup for Navier–Stokes**.

> **Goal:** reproduce the paper's velocity construction as an auditable computational object, from the similarity coordinates and leading axisymmetric vortex through the all-order background corrections, oscillatory stress realization, local field, and final compact localization.

## What is implemented now

The repository starts with the exact kinematic layer used by the paper:

- similarity coordinates from Eq. (4.1), including numerical solution of the implicit concentration scale `q(z,t)`;
- the leading axisymmetric field from Eqs. (4.3)–(4.7);
- exact cylindrical-to-Cartesian conversion corresponding to Eq. (4.5);
- the formal all-order background expansion pattern from Eq. (5.1);
- explicit data structures for the summed local representation `u_loc = curl(A) + B e_theta` from Eq. (9.21) and localization `u = curl(c A) + c B e_theta` from Eq. (10.4);
- numerical sanity checks for coordinate identities, incompressibility, and the predicted blow-up scaling.

The difficult profile and correction constructors are deliberately **not replaced by guessed closed forms**. They are tracked as reconstruction stages in `docs/RECONSTRUCTION_PLAN.md`. A result is only labeled `paper-exact` after its choices and identities are tied to the paper/official Lean formalization.

## Core equations

Let

```text
tau = 1 - t
A   = 1/2 + h
D   = 1/2 - h
z   = q^D eta
tau = q (1 - eta^2)
X   = r^2 / (2q)
```

with `0 < h < 1/100`. Equivalently, `q` is the unique positive solution of

```text
q - z^2 q^(2h) = tau.
```

For leading profiles `E(X,eta)` and `U(X,eta)`, define

```text
A_X(U) = (1/X) integral_0^X U(x,eta) dx
d       = 1 - eta^2
L       = 1 - 2 h eta^2
V0      = X/L * [2 eta U - 2 D eta A_X(U) - d d_eta A_X(U)]
```

and

```text
u_theta^(0) = q^(-A) E
u_z^(0)     = q^(-A) U
r u_r^(0)   = V0.
```

At `z=0`, `q=tau`, and along `r=sqrt(2 X_in tau)` the paper proves

```text
u_theta = tau^(-A) (e0 + O(tau^(2h))) -> +infinity.
```

## Install

```bash
python -m pip install -e .
python -m pytest
```

Run the leading-field demo:

```bash
python examples/leading_field_demo.py
```

## Source of truth

- OpenAI paper: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- OpenAI Lean formalization: https://github.com/openai/NavierStokesAndEuler
- OpenAI announcement: https://openai.com/index/navier-stokes-solution/

This repository is an independent reconstruction and is not an OpenAI repository.
