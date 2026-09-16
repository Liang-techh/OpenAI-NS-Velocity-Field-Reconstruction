# First-Picard wide pressure source and chain provenance

Status: **formal-structure only**.

This increment binds the pinned pressure source and pressure chain to the
genuine theorem-selected first Picard state `x1`:

```text
source = a^2 * (phi1 * phi1)

pressure = j1(-4*A*eta*primitive(source)
              + d*parameterPrimitive(source)
              - 2*eta*mulY(source)).
```

The existing first-Picard angular-square state supplies five ordinary channels
for `phi1 * phi1`, with powers `Lambda^0` through `Lambda^-4`.  The new state
returns normalized Decimal five-tuples `C_0,...,C_4` satisfying

```text
d_eta^m source_n(eta) = a(eta)^2 * sum(C_p[n,m](eta)/Lambda^p),
d_eta^m pressure_n(eta) = a(eta)^2 * sum(P_p[n,m](eta)/Lambda^p).
```

The amplitude-square Bell factors are reused from the landed wide natural
source.  Thus positive eta derivatives include the derivatives of the actual
theorem-selected amplitude; the implementation never differentiates a
zeroth-order normalized value while ignoring those terms.  The common `a^2`
scale stays outside the Decimal tuples, and `source_jet_log` /
`pressure_jet_log` expose signed-log views when the full magnitude is needed.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main files:
  - `NavierStokes/AxisContraction.lean`
  - `NavierStokes/AxisOperators.lean`
  - `NavierStokes/AxisWeightEstimates.lean`
  - `NavierStokes/NaturalAxisCoefficients.lean`
- Relevant symbols:
  - `NavierStokes.AxisContraction.naturalRemainder`
  - `NavierStokes.AxisOperators.productFamily`
  - `NavierStokes.AxisOperators.primitive`
  - `NavierStokes.AxisOperators.parameterPrimitive`
  - `NavierStokes.AxisOperators.mulY`
  - `NavierStokes.AxisOperators.regularInverse`

The exact upstream `naturalRemainder` source and pressure expression are at
`AxisContraction.lean:339-365`.  The coefficient product and eta-Leibniz law
are at `AxisOperators.lean:86-95`; the radial primitive, parameter primitive,
and multiplication-by-`Y` row maps are supplied by the operator bindings at
`AxisOperators.lean:681-717`.  The regular inverse divisor is pinned by
`AxisWeightEstimates.lean:388`, namely
`radialDivisor(1,n-1) = n^2` for output row `n >= 1`.  The coefficient window
is `[-11/10,11/10]` from `NaturalAxisCoefficients.lean:27`.

The executable prerequisites are the actual first-Picard angular square in
`src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_angular_square.py`,
the x0 wide source in
`src/openai_ns_reconstruction/axis_coefficient_wide_natural_source.py`, and
the x0 pressure row maps in
`src/openai_ns_reconstruction/axis_coefficient_wide_natural_pressure.py`.
The implementation itself is
`src/openai_ns_reconstruction/axis_coefficient_wide_first_picard_pressure.py`.

## Executable binding and scale boundary

`wide_first_picard_pressure_state(x1)` accepts only an
`ActualScheduleWideFirstPicardState`.  It rebuilds fixed `AxisData` from
`x1.reference`, reuses the amplitude and wide-source chain carried by
`x1.remainder.axial.wide_pressure`, and checks the shared `epsilon` and
theorem-selected `Lambda`.  The public methods
`normalized_source_factors(n,m,eta)` and
`normalized_pressure_factors(n,m,eta)` return five Decimal numerators in
power order `p=0,...,4`.  Their aliases `source_jet` and `pressure_jet` expose
the same normalized tuples; the `*_jet_log` methods attach the common
`2*log(a(eta))` scale and the explicit `p*log(Lambda)` factor.

All radial sums are finite and exact at the requested row.  Eta products use
the finite binomial Leibniz sum.  `primitive` reads the predecessor row and
divides by `n`, `parameterPrimitive` reads the predecessor at eta order
`m+1` and divides by `n`, `mulY` reads the predecessor without division, and
the final `j1` reads the pressure input at row `n-1` and divides by `n^2`.
The row-zero pressure output is exactly zero after the pinned regular inverse.
No outer `inverseL` multiplication is included in this raw pressure chain.

## Truth boundary

The state exposes the actual five-channel x1 pressure source and pressure chain
while retaining `paper_exact = false`,
`natural_remainder_x1_materialized = false`, `fixed_point_materialized =
false`, and `global_axis_norm_certified = false`.  This does not materialize
the full `naturalRemainder(x1)`, the next Picard iterate, a fixed point, the
derived final average/pressure fields, `NaturalProfileAssembly`, global
weighted `AxisSpace` bounds, or the paper-exact velocity profile.  The actual
phase and amplitude values also retain the landed numerical quadrature
boundary.
