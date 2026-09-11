# Uniform LocalBase → rounded-normal provenance

Status: **formal-structure**. This note does not certify a paper-exact base field or a completed oscillatory correction.

## Source mapping

Pinned upstream source: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, primarily `NavierStokes/PhaseEstimates.lean`.

The executable adapter in `src/openai_ns_reconstruction/phase_uniform_bounds.py` mirrors the quantitative part of `PhaseEstimates.LocalBaseBounds` and the theorem `localBase_derivative_errors` used immediately before `rounded_normal_estimates`.

For a convex slow-domain `U`, actual fields `F,G`, order-zero fields `F0,G0`, and a uniform constant `M`, the pinned structure contains, among its hypotheses:

- second derivative bounds on `F0,G0` by `M`;
- a first derivative bound on `F0` by `M`;
- `|F-F0| <= M epsilon^2`;
- first-derivative actual/reference errors for `F,G` bounded by `M epsilon^2`;
- actual directional bounds `|slowR F|, |slowZ F|, |slowZ G| <= M`.

The pinned theorem `localBase_derivative_errors` then gives, for `q,q0 in U` and `||q-q0|| <= d`,

`|slowR F(q)-slowR F0(q0)| <= M (d + epsilon^2)`

and the identical bound for `G`. Hence a separately certified slow-box diameter `d <= S^-3` implies exactly the two inputs

`|FR-FR0|, |GR-GR0| <= M (S^-3 + epsilon^2)`

required by the already-landed `PhaseEstimates.rounded_normal_estimates` adapter. The LocalBase directional hypotheses also directly provide the uniform `|FR|, |FZ|, |GZ| <= M` inputs.

## What the executable module certifies

`LocalBaseNumericEnvelope` is a fail-closed holder for **already certified uniform sup bounds**. It rejects entries above the exact numeric thresholds in the pinned LocalBase structure. It intentionally does not accept samples or estimate suprema.

`UniformLocalBaseBridge` requires exact agreement with the existing `PhaseScaleCertificate` in `M` and `epsilon`, and requires a certified domain diameter `d <= S^-3`. It exposes both the sharper `M(d+epsilon^2)` error and the theorem-facing `M(S^-3+epsilon^2)` envelope.

Focused tests independently instantiate simple polynomial derivative families to cross-check the mean-value arithmetic. Those polynomials are test fixtures only and are not claimed to be any field from the paper.

## Boundary / remaining obligations

This increment does **not** materialize Proposition 5.5's paper-exact base functions `F,G,F0,G0`; it does not infer convexity, differentiability, second derivative bounds, or actual/reference errors from numerical samples. Those hypotheses must come from the still-unfinished leading-profile/all-order-background construction or an independent analytic/formal certificate.

It also does not by itself prove all of Eqs. (7.9)-(7.11). It only closes the LocalBase/C1-C2 implication needed to upgrade the existing pointwise rounded-normal adapter once genuine base-field bounds are supplied. Representative-frequency/inverse-radius/slot hypotheses, moving-frame and damping estimates, stress-cone decomposition, amplitudes, curl waves, mean corrections, and Section 9 residual iteration remain separate obligations.

Accordingly Stage 3 remains `formal-structure`, and `paper_exact_velocity_available` must remain `false`.
