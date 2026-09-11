# Analytic-input norm certificate provenance

Status: **formal-structure only**. This file does not certify a paper-exact leading profile and does not change `paper_exact_velocity_available`.

## Pinned source mapping

Pinned official source commit: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The executable adapter in `src/openai_ns_reconstruction/axis_analytic_input_bounds.py` follows these exact Lean constructions:

- `NavierStokes/NaturalAxisCoefficients.lean`
  - `exists_coefficientFamily_on_neighborhood`: the constructed coefficient radius is `epsilon = rho / 2` and the common coefficient-space bound is `B * radiusLoss ((rho/2)/rho)`.
  - `exists_coefficientFamily_with_convex_domain` / `exists_analyticInputs`: the same `rho` is retained as the larger analytic radius, so the ratio is exactly `1/2`.
  - `AnalyticInputs.amplitudeBound`: `radiusLoss (coefficients.epsilon / radius)`.
  - `CoefficientFamily.norm_le`: every fixed field used by `CoefficientFamily.axisData` has norm bounded by the same family bound.
- `NavierStokes/AnalyticCoefficientBounds.lean`
  - `radiusLoss(q) = sum_{m>=0} (m+1)^2 q^m`. For `q=1/2`, the elementary closed form `(1+q)/(1-q)^3` gives the exact value `12`.
- `NavierStokes/AxisResolvent.lean`
  - `factorialMajorant K k = K^k / (k! (k+1)!)`.
  - `naturalOperator_pow_bound` uses `K = 2560 * ||chi||`.
  - `naturalResolvent_norm_le` bounds the resolvent norm by `sum_k factorialMajorant (2560*||chi||) k`.

Therefore, once a positive common analytic radius `rho` and a certified common complex-field value bound `B` are supplied for the actual compact neighborhood chosen in the coefficient-family construction, the adapter derives:

- `epsilon = rho/2`;
- common fixed-field norm upper `12 B`;
- normalized angular-amplitude norm upper `M = 12`;
- `||chi|| <= 12 B`;
- resolvent majorant parameter `K <= 30720 B`;
- an executable upper enclosure for the complete factorial resolvent series.

The factorial-series evaluator uses the exact positive-term recurrence `a_(k+1)/a_k = K/((k+1)(k+2))`. After this ratio is at most `1/2`, the untouched tail is bounded by a geometric series. Decimal arithmetic is rounded toward `+infinity`, followed by an upward binary64 conversion. Tests independently cross-check the result against the Bessel identity `sum K^k/(k!(k+1)!) = I_1(2 sqrt(K))/sqrt(K)`.

## What this closes

This removes two previously opaque caller choices from the landed `axis_remainder_bounds.py` path for the canonical Lean coefficient-family construction: the normalized amplitude bound is no longer arbitrary (`M=12`), and the natural-resolvent norm can be generated from the single certified common complex-field bound `B` rather than being supplied as an unrelated constant. The fixed `AxisData` coefficient-space norm ledger also collapses to the common `12 B` upper bound.

## What remains open

This increment deliberately does **not** manufacture the missing analytic witnesses. In particular it does not yet construct an explicit admissible `rho` for the actual schedule pressure, does not derive the true common complex-field sup bound `B` on the corresponding compact set, and does not certify `realPartSup(axisPhase, compactSet)`. It therefore does not yet instantiate the numerical fixed-point contraction or materialize `phi/u/average/pressure`.

The resulting bounds can also become too large for downstream binary64 remainder propagation; conversion fails closed when a theorem-valid real bound exceeds binary64 range. A later step may need arbitrary-precision/log-domain propagation, but this file does not weaken or replace the theorem inequalities to avoid that issue.
