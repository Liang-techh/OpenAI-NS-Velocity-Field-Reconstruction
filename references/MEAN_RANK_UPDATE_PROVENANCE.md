# Section 8 five-row mean-rank update provenance

Status: **formal-structure**. `paper_exact = false`.

## Pinned source

This increment is mapped to the paper's Section 8 / five-equation mean repair
(Equation (35) in the pinned formalization) and to
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/FiveRowRank.lean`
  - `Debt = (P, Jθ, Jz)`;
  - `cellStep`, `cellLower`, `cellUpper`;
  - angular powers `(2, -2-2λ, -2λ)`;
  - axial powers `(1, 1-2λ)`;
  - normalized angular debt `(0, -d_P/(2C), d_Jz/C)`;
  - normalized axial debt `(0, -d_Jθ/C)`;
  - compact `deltaV` / `gamma` moment repair and the five-row theorem;
- `NavierStokes/LocalizedMomentRepair.lean`
  - the support of each constructed bump lies in
    `((3l+u)/4, (l+3u)/4)` inside its separated cell;
- `NavierStokes/MeanRankUpdate.lean`
  - `scaleField`, `scaleDebt`, `normalizeDebt`;
  - the physical angular/axial increments;
  - the physical five-row scaling identities and support inside `(ell*a,ell*b)`.

## Executable artifact

`src/openai_ns_reconstruction/mean_rank_update.py` now exposes the pinned
finite-dimensional structure without requiring the unresolved Section 6--7
paper field:

1. it constructs the exact separated three-cell angular and two-cell axial
   geometry and the exact inner support intervals used by the Lean repair;
2. it implements the official physical debt normalization
   `(d_P/U^2, d_Jθ/(ell^3 U^2), d_Jz/(ell^2 U^2))`;
3. it builds the angular/axial power-moment systems and the corresponding
   target moments;
4. it applies the physical scaling identities, so a normalized five-moment
   repair maps to rows `(0,0,-d_P,-d_Jθ,-d_Jz)`;
5. for executable evaluation it solves those two finite moment systems and
   returns compact angular/desired-axial increments on the correctly scaled
   support intervals.

The executable bump used for step 5 is the repository's explicit C-infinity
`CompactMomentBump`, restricted to the same official inner support intervals.
The pinned Lean witness instead uses Mathlib's noncomputable `ContDiffBump`.
Their transition/interior values are **not** claimed to coincide. Numerical
Gauss--Legendre moments and float64 matrix inversion are therefore execution
infrastructure, not a formal proof of the five-row theorem.

## Independent regression

`tests/test_mean_rank_update.py` checks the cell geometry and debt powers by
hand and then verifies the completed physical five rows by a separate SciPy
adaptive quadrature over the physical support intervals. The verification does
not reuse the production moment-matrix integrator. The test data are synthetic
fixtures used only to cross-check the algebra; they are not manuscript wave or
background data.

## Remaining paper-exact boundary

This increment does **not** supply the actual mean debt, `λ,C,ell,U`, or the
power-law patch from the completed Sections 6--8 construction. Those values
must be derived from the genuine background, pulse/covariance and oscillatory
wave hierarchy. It also does not prove that the executable bump may replace
the pinned noncomputable Mathlib bump, nor provide interval/formal certificates
for its moment matrix and integrals.

Consequently this module is a reusable compact-mean solver/validator only. It
does not construct a surrogate wave, does not close Section 8, does not start
the Section 9 residual-improvement iteration, and must not change
`paper_exact_velocity_available` from `false`.
