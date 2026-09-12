# Stage-1 mixed-scale axial natural remainder provenance

Status: **formal-structure**. This increment does not claim a paper-exact leading profile or a materialized fixed point.

## Pinned source

Mapped to `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, principally:

- `NavierStokes/AxisContraction.lean`
  - `naturalRemainder`
  - the axial `lin2`, `slow2`, `pressure` branch
- `NavierStokes/AxisOperators.lean`
- `NavierStokes/NaturalAxisCoefficients.lean`

At the contraction centre `x0 = referencePair`, the pinned axial remainder has the form

`inverseL * (lin2 - (1/Lambda) * slow2 + pressure)`.

The pressure source is the theorem-selected `a^2 * phi0^2`, with
`a = exp(Lambda * realPhase) / C`.

## Landed capability

`src/openai_ns_reconstruction/axis_coefficient_wide_axial_remainder.py` binds this formula to the actual `SchedulePressure` reference state and the already-landed theorem-selected wide amplitude/pressure chain.

The module intentionally preserves three scales instead of collapsing them:

1. `inverseL * lin2` remains an ordinary coefficient state;
2. `-inverseL * slow2 / Lambda` stores the ordinary coefficient as an exact `Decimal.from_float` numerator together with the actual theorem-selected Decimal `Lambda` denominator;
3. `inverseL * pressure` keeps the genuine nonzero `a^2` magnitude in `SignedLogCoefficientJet` form while applying the last `inverseL` multiplication with the exact finite radial convolution and eta-Leibniz rule.

There is no caller override for pressure data, sigma, epsilon, `Lambda`, `C`, amplitude, source coefficients, or the coefficient operators. No radial cutoff or fitted derivative table is introduced.

## Independent verification

The regression suite cross-checks the ordinary pieces against the previously landed full binary64 `naturalRemainder` composition using regression-only isolating values:

- `a=0, t=0` isolates `inverseL * lin2`;
- `a=0`, the difference between `t=0` and `t=1` isolates `inverseL * slow2`.

Those isolating values are test oracles only and are never promoted to manuscript parameters.

For the pressure branch, eta derivative order zero is independently compared with the existing full remainder evaluated at the exact constant coefficient state `a=1`; subtracting the `a=0` result isolates the complete `inverseL * pressure` chain. Higher eta derivatives remain represented through the already-verified wide pressure Bell/Leibniz recurrence plus the new exact `inverseL` Leibniz convolution.

The tests also require a nonzero pressure coefficient to remain nonzero in signed-log form and require binary64 projection to fail closed rather than silently underflow to zero.

## Explicit non-claims / next blocker

This is **not yet the complete `naturalRemainder(x0)`**. The angular branch still needs the same theorem-selected `1/Lambda` scale treatment through `naturalResolvent`, and the outer Picard factor `1/(2*Lambda)` has not yet been applied to the mixed-scale pair. Consequently no genuine `x1`, fixed-point `phi/u`, radial average/pressure reconstruction, `NaturalProfileAssembly`, global all-index `AxisSpace` certificate, or support/moment/matching/cone closure is claimed.

`paper_exact_velocity_available=false` remains mandatory.
