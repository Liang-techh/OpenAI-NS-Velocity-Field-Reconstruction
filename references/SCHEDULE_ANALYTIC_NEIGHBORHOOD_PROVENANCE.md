# Schedule analytic-neighborhood provenance

Status: **formal-structure**. This file records an executable conservative realization of the compactness choices used by the pinned natural-axis coefficient construction. It does **not** make Stage 1 paper-exact and does not construct the coefficient-space fixed point.

Pinned official source commit: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

## Official source mapping

- `NavierStokes/PressureDatum.lean`
  - `strip := {z | |z.im| < 1/2}`.
  - `complexPressure_analytic` proves the pressure extension is analytic on that strip.
- `NavierStokes/SchedulePressure.lean`
  - `shapeExponent_bounds` gives `0 <= shapeExponent <= 1`.
  - `admissible` instantiates `PressureDatum.Admissible ... 1` for the actual `clockWeight` / `shapeExponent` schedule.
  - `complexAxisPressure_analytic` identifies the actual schedule pressure with the above holomorphic extension.
- `NavierStokes/NaturalAxisCoefficients.lean`
  - `regularSet` is the intersection of the pressure strip, `complexL != 0`, and `denominator != 0`.
  - `exists_common_neighborhood` chooses a positive common radius nonconstructively by compactness.
  - `finite_family_bound` chooses one positive common complex-field bound in the form `1 + sum_i |b_i|`.
  - `exists_coefficientFamily_on_neighborhood` then sets the coefficient-space radius to `epsilon = rho/2`.
  - `axisPhase` is the primitive of `complexGradient`.

## Executable replacement of the opaque compactness choices

`src/openai_ns_reconstruction/schedule_analytic_neighborhood.py` makes three previously opaque upstream quantities explicit for the already-landed actual schedule.

### 1. Common complex radius `rho`

The real parameter window is the pinned `[-11/10,11/10]`. We first restrict to radius `1/4`, exactly half the available pressure-strip half-width, so every point in the resulting tube satisfies `|Im z| <= 1/4 < 1/2`.

For `complexL(z)=1-2 h z^2`, the real window has the explicit positive lower bound

`L_min = 1 - 2 h (11/10)^2`.

A derivative bound on the fixed search tube keeps `complexL` away from zero.

For the second rational denominator, use the exact factorization

`H(z)^2 + sigma^2 = (H(z)-i sigma)(H(z)+i sigma)`.

At every real window point `x`, `H(x)` is real, hence each factor has modulus at least `sigma`. A polynomial derivative upper bound `M_H'` on the fixed search tube gives

`|H(z)-H(x)| <= M_H' |z-x|`.

The selected radius enforces `M_H' rho <= sigma/4`; therefore both factors stay at least `3 sigma/4`, giving a strictly positive denominator throughout the certified tube. This avoids the much weaker `O(sigma^2)` radius that would result from perturbing the product as a whole.

The actual schedule sigma is not caller-chosen: `certify_actual_schedule_analytic_inputs` reuses the theorem-side sigma already produced by `schedule_axis_margin.certify_schedule_low_Z_margin`.

### 2. Actual finite-family complex bound `B`

On the certified tube, `|1+z^2| = |z-i||z+i| >= (1-rho)^2`. Since the actual schedule has exponent cap one, the complex pressure kernel obeys

`|complexKernel(a,z)| <= (1-rho)^(-4)`.

Combining this with the existing no-sampling `clock_mass_upper` gives explicit bounds for the actual complex pressure and its first derivative. Those bounds, together with the polynomial/rational bounds above, give one bound for each of the eleven pinned `NaturalAxisCoefficients.Field` values:

`one, eta, d, inverseL, uStar, uStarEta, wStar, hStar, zStar, chi, gradient`.

For `chi` and `gradient` the factor separation is used directly. If

`c = sup |H(z)-H(x)| / sigma <= 1/4`, then

`|H^2/(H^2+sigma^2)| <= 1 + (1-c)^(-2)`

and

`|H/(H^2+sigma^2)| <= 1 / (sigma (1-c))`.

This is substantially tighter than bounding the numerator and denominator independently.

The common family bound is then chosen in the same constructive shape as pinned `finite_family_bound`:

`B = 1 + sum(field_bounds)`.

No grid maximum or fitted complex-field value is used to construct `B`.

### 3. Complex phase real-part supremum

The certified closed tube is convex and contains zero. Since `axisPhase` is a primitive of `complexGradient`, the straight-line integral from zero to any `z` in the tube yields

`|axisPhase(z)| <= |z| sup_tube |complexGradient|`.

The module records the resulting finite upper bound as `axis_phase_real_part_sup_upper`, which is also a valid upper bound for `sup Re(axisPhase)` required by the landed `natural_scale_selection.py` normalization-threshold adapter.

## Independent regression boundary

`tests/test_schedule_analytic_neighborhood.py` samples complex points only as a consequence check. The radius and field bounds are built analytically before any sample is taken. The tests directly evaluate the polynomial/rational `complexL`, `H`, `chi`, and `gradient` formulas inside the certified tube, and independently cross-check the pressure / pressure-derivative / `zStar` consequences against the existing actual schedule evaluator on the real axis. Real-phase samples are also checked against the analytic phase-sup envelope.

## Remaining Stage-1 blockers

This increment removes the need for a caller-supplied common neighborhood radius, common complex-field sup bound, and compact-set phase real-part bound for the actual schedule. It does **not** yet:

- materialize the coefficient-space fixed-point fields `phi/u/average/pressure`;
- prove the Python real-arithmetic inequalities as interval or Lean proof objects;
- address the potentially enormous coefficient-space norm / normalization scales induced by the conservative theorem-side sigma;
- connect a solved fixed point to `NaturalProfileAssembly`;
- verify the final support, moments, matching, and cone conditions.

Accordingly `paper_exact_velocity_available` must remain `false`.
