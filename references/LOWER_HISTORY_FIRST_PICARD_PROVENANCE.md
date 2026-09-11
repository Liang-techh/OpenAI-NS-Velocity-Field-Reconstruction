# Section 5 lower-history to first PositiveAxis Picard step provenance

Status: **formal-structure**.

## Pinned official source

This increment is mapped to
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd` and the paper's Section 5, especially Eqs. (5.3)-(5.8):

- `NavierStokes/PositiveAxisSystem.lean`
  - `actualLowerSource`
  - `baseAtOrderZero`
  - the displayed `A0`, `A1`, and forcing definitions;
- `NavierStokes/PositiveAxisExistence.lean`
  - `lowerHistoryData`
  - the actual symmetric Volterra solution / singular first-order system bridge;
- `NavierStokes/SlowRecursion.lean`
  - the induction step supplying only already-constructed orders `< n`;
- the manuscript Lemma 5.1 / Eq. (5.7), where
  `W_n = sum_{k>=0} K^k G f_n` and the first term is `G f_n`.

The Eq. (5.7) radial variable is `xi=sqrt(X)`.  Therefore a lower coefficient
jet used at a Volterra quadrature point `xi=s` must be evaluated at the squared
radius `X=s^2`; replacing this by an interpolation in `xi` would not be the
pinned coefficient system.

## Executable artifact

`src/openai_ns_reconstruction/background_lower_history_solver.py` connects the
already-landed exact layers instead of introducing a new surrogate solve:

1. `lower_history_positive_axis_providers` queries only orders
   `0,...,n-1` and maps every Eq. (5.7) radius `xi` to `X=xi^2` before calling
   `positive_axis_point_data_from_lower_history`;
2. `positive_axis_eq_5_7_fields_from_lower_history` feeds those derived
   `BaseJet` / `SourceJet` values into the pinned `A0/A1/f_n` formulas;
3. `first_positive_order_term_from_lower_history` evaluates the genuine first
   Lemma-5.1 term `G f_n` over the whole radial interval, so the forcing at the
   quadrature nodes is rebuilt from strict lower history, preceding double-`Z`
   diffusion, and the regular Eq. (5.6) `Omega/X` source.  It no longer accepts
   an unrelated six-vector forcing callback at this boundary.

## Independent check

`tests/test_background_lower_history_solver.py` uses an order-one polynomial
lower history that is independent of `eta` and evaluates the production
`G f_1` bridge.  The oracle does **not** call the production lower-source or
singular-inverse routines: at `eta=0` it hand-reduces the two composed `Z`
operators and Eq. (5.6) to explicit polynomial formulas, then integrates each
weighted component independently with SciPy adaptive quadrature.  A separate
regression records every provider query and verifies that only strict lower
orders are requested and that the radial mapping is exactly `X=xi^2`.

These polynomial data are verification fixtures only; they are not claimed to
be the paper's coefficient hierarchy.

## Boundary / remaining blocker

This bridge performs one real theorem-shaped positive-order solver step
**conditional on** lower-history second-jet functions.  It does not construct
those lower jets.  For paper-derived recursion they must come from the genuine
Section-4 leading profile and, at positive lower orders, from the converged
Lemma-5.1 solution after Lemma-5.2 compact repair.  Sampled or fitted jets are
not a substitute.

The module also does not certify the common complex strip or the actual
Picard-tail constants `C_n`, does not sum the full Picard series, does not derive
the hierarchy's moment/patch-factor jets, does not prove the infinite cutoff
schedule/local finiteness, and does not establish Proposition 5.3 all-jets-flat
residual decay.

Therefore Stage 2 remains **formal-structure** and
`paper_exact_velocity_available` remains `false`.
