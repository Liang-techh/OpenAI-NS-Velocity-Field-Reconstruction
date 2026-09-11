# Positive-order axis system provenance

Status: **formal-structure**.

## Pinned source

This increment is a direct executable transcription of the real-valued algebra in
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`,
`NavierStokes/PositiveAxisSystem.lean`.

The pinned file defines `Jet`, `BaseJet`, `SourceJet`, the six singular weights
`(0,0,2,0,3,1)`, the displayed matrices `A0` and `A1`, and the forcing vector used
by the positive-order system. Its theorem `matrixRHS_jet` expands direct matrix
multiplication back into the four profile equations, including the pressure-row
substitution. `A1_shape` and `A1_high_parameters_irrelevant` record the sparse
parameter-derivative structure used in the Lemma 5.1 derivative count.

## Executable artifact

`src/openai_ns_reconstruction/background_positive_axis.py` implements:

- structured `RadialParameterJet`, `PositiveAxisBaseJet`, and
  `PositiveAxisSourceJet` inputs;
- the exact real `A0(h, lambda_n, C, xi, eta, base)` formula at
  `lambda_n = 2 n h`;
- the exact sparse `A1(h, xi, eta, base)` formula;
- the exact forcing vector
  `(0,0,0,2 xi pressureSource,2 source.angular,
   2 source.axial - 4 eta xi^2 pressureSource / ell)`;
- a typed `positive_axis_eq_5_7_fields` adapter into the already-landed
  Eq. (5.7) Picard solver; and
- `first_positive_order_term_from_jets`, which evaluates the genuine first
  Lemma-5.1 term `G f_n` after constructing `f_n` from a structured lower-order
  `SourceJet`, rather than accepting an arbitrary six-vector forcing callback.

`tests/test_background_positive_axis.py` independently evaluates the expanded
profile right-hand sides and compares them with `A0 W + A1 d_eta W + f_n`.
The test does not call production scalar helpers for that oracle. A second test
uses a constant `SourceJet`, for which the three nonzero components of `G f_n`
have closed-form polynomial integrals, and compares those formulas against the
actual singular inverse. The sparse `A1` shape and fail-closed validation are
covered separately.

## What this closes

The repository no longer needs arbitrary caller-supplied **matrices** in order
to exercise the Eq. (5.7) solver. Given theorem-shaped base and lower-history
jets, the exact matrix/source algebra and a true positive-order `G f_n` step are
now executable. This is the concrete Eqs. (5.3)-(5.6) -> Eq. (5.7) coefficient
interface that was previously missing between `background_recurrence.py` and
`background_inner_solver.py`.

## Boundary / blocker

This does **not** make Stage 2 paper-exact. The `BaseJet` and `SourceJet`
providers are still inputs. In the paper/Lean recursion they must be populated
from the materialized Section-4 leading profile and the already-finalized lower
orders (including the regular previous `Omega/X` source). Issue #1 has not yet
materialized the complete leading fixed-point profile, and this increment does
not manufacture lower-history convolutions from sampled or toy data.

It also does not certify the common complex strip, the constants `C_n` used in
Eq. (5.8), convergence of the full Picard series for the actual hierarchy, the
Lemma 5.2 uniform nonvanishing/repair hypotheses, or the infinite cutoff/residual
flatness conclusions. Numerical evaluation is binary64 and the radial inverse
uses the repository's finite Gauss-Legendre rule.

Therefore Stage 2 remains **formal-structure** and
`paper_exact_velocity_available` must remain `false`.
