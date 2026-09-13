# Generic mixed-scale coefficient algebra provenance

Status: **formal-structure only**.

Artifact:

- `src/openai_ns_reconstruction/axis_coefficient_mixed_scale.py`
- `references/provenance_manifest_addendum_axis_coefficient_mixed_scale.json`

This increment supplies the reusable sparse channel algebra needed by later
Picard layers.  A `MixedScaleCoefficient` is an immutable map
`(q, p) -> Decimal`, representing the formal channel
`a(eta)^q Lambda^-p` by its numerator alone.  The implementation accepts
nonnegative integer amplitude powers `q` (including odd `q` for the physical
ratio `F = a phi`) and nonnegative inverse-Lambda powers `p`; zero numerators
are omitted from the sparse support, and every nonzero channel is retained
without a fixed-width truncation.

The stored numerator has the full eta derivative convention: it is the actual
eta derivative of the complete channel divided by the common `a(eta)^q`
factor.  The family eta-derivative primitive therefore requests provider order
`m + 1`.  It does not differentiate a stored numerator in place.  A bare
amplitude-power shift is intentionally absent.  To multiply by an amplitude
factor, callers must provide its full Bell derivative family and use
`mixed_scale_product`, which applies the eta Leibniz rule to both families.
Lambda shifts are available because Lambda is eta-independent.

The public family primitives are:

- `mixed_scale_product(left, right)`, using radial convolution over `i + j = n`
  and the binomial eta sum over `k + l = m`;
- `mixed_scale_eta_derivative(source)`, using the source jet at `m + 1`;
- `mixed_scale_euler(source)`, multiplying a row by `n`;
- `mixed_scale_average(source)`, dividing a row by `n + 1`;
- `mixed_scale_primitive(source)`, returning zero at `n = 0` and
  `source(n - 1, m) / n` otherwise;
- `mixed_scale_j_r(source, r)`, for `r = 1` or `2`, returning zero at `n = 0`
  and `source(n - 1, m) / (n * (n + r - 1))` otherwise.

The zero rows of the primitive and `J_r` return structurally without probing
the source row.  This preserves triangular causality for providers whose
radial recursion depends on the output row at `n = 0`.

All arithmetic runs inside a local 96-digit Decimal context.  No amplitude or
Lambda power is numerically formed, and no binary64 scalar conversion is
accepted by the channel operations.

## Pinned formal source

The algebra follows the pinned `openai/NavierStokesAndEuler` source at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/AxisOperators.lean:86-95` for the coefficient product;
- `NavierStokes/AxisOperators.lean:501-504` for the regular inverse family;
- `NavierStokes/AxisOperators.lean:512-518` for the differential-family
  parameter/Euler actions;
- `NavierStokes/AxisOperators.lean:681-717` for the coefficient-operator
  bindings;
- `NavierStokes/AxisWeightEstimates.lean:388` for the radial divisor
  `n * (n + r - 1)`;
- `NavierStokes/NaturalAxisCoefficients.lean:27` for the coefficient-domain
  window inherited by actual providers.

The current repository base remains main commit
`67ab5ecd32f80225df8bce2760b943c20a7d3d7b`.  This module is independent of
the reviewed PR 279 additions and does not alter them.

## Truth and dependency boundary

This is a generic typed arithmetic layer.  It does not select actual
SchedulePressure data, provide an arbitrary Picard evaluator, materialize
`naturalRemainder(x1)` or `naturalRemainder(x2)`, choose an amplitude or
Lambda, establish a global weighted `AxisSpace` bound, or certify convergence
to a fixed point.  Passing a hand-built provider through these functions is an
algebra regression only; it does not promote that provider to paper data.

The provenance status therefore remains `formal-structure`, with
`full_reconstruction=false` and `paper_exact_velocity_available=false`.
