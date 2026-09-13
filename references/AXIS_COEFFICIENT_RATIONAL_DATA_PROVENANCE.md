# Rational `NaturalAxisData` input provenance

Status: **exact finite scalar algebra only**.  This artifact does not promote
the coefficient family, the pressure datum, or the velocity to paper-exact.

## Pin and scope

The official pin is `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.  The definitions are in
`NavierStokes/NaturalAxisData.lean:24-35`; the field enumeration and real
gradient are in `NavierStokes/NaturalAxisCoefficients.lean:23-49,105-134`;
the angular recurrence and closed form are in
`NavierStokes/AxisReference.lean:61-98`.

For chosen binary64 values `h`, `j`, `sigma`, and `eta`, interpret each input
as the exact dyadic `Fraction.from_float(value)`, rather than as its displayed
decimal.  Require the theorem-side signs and window (`0 < h,j <= 1/1000`,
`sigma > 0`, and `eta in [-11/10,11/10]`); reject a zero `L` or denominator.

## Exact fields and Taylor rule

The scalar fields are exactly `A = 1/2 + h` and `D = 1/2 - h`.  For a formal
shift `t`, set

```text
e = eta + t                 d = 1 - e^2
U = 4*e + j                 L = 1 - 2*h*e^2
H = D*e + d*U               W = 1 - 4*d - 2*D*e*U
G = H^2 + sigma^2
chi = H^2/G                 normalizedGradient = -L*H/G.
```

The nine fixed fields are `one = 1`, `eta = e`, `d`, `inverseL = 1/L`,
`uStar = U`, `uStarEta = 4`, `wStar = W`, `hStar = H`, and
`normalizedGradient`; all have radial rows `n > 0` equal to zero.  Product
Taylor coefficients use the exact Cauchy rule
`(fg)[r] = sum(k=0..r) f[k]*g[r-k]`.  Each reciprocal uses
`b[0]=1/a[0]` and `b[r]=-sum(k=1..r) a[k]*b[r-k]/a[0]`.

For every finite order, build `G^-1` by that recurrence, then form `chi` and
`normalizedGradient` by exact Fraction products.  A Taylor coefficient `q[m]`
represents the actual eta derivative `m! * q[m]`.

The pinned angular reference is independently recovered with
`c = -chi/2`, `phi0[0] = 1`, and
`phi0[n+1] = c*phi0[n]/((n+1)*(n+2))`; equivalently

```text
phi0[n] = (-chi/2)^n / (n! * (n+1)!).
```

Thus all finite `(n,m,eta)` angular reference derivatives are exact rational
outputs of the Fraction Taylor algebra, subject only to the stated nonzero
denominators.

## Existing implementation comparison

The cached local definitions are `natural_axis.py` (`D`, `A`, `d`, `L`,
`axis_U`, `H`, `W`, `chi`, `real_gradient`),
`axis_coefficient_data.py` (`_base_taylors`, `_field_taylor`), and
`axis_reference_angular_jets.py` (`_h_taylor`, `_chi_taylor`,
`normalized_taylor`).  No algebraic formula discrepancy was found: the old
implementations are binary64 evaluations of the same expressions.  They can
differ in last bits because rounding and, for `chi`, cancellation handling occur
during float operations; they are point values, not exact enclosures.

Converting each exact Fraction result `q` to a directed Decimal interval at
precision `p` means `[ROUND_FLOOR_p(q), ROUND_CEILING_p(q)]`.  This enclosure
applies only to the rational scalar fields and finite angular-reference Taylor
layer.  It does not enclose the old binary64 computation, accumulated
coefficient arithmetic, the `zStar` pressure chain, or the amplitude/phase
oracle.

## Boundary

`zStar = Z(h,j,P)` still depends on the actual SchedulePressure datum and its
derivative.  `realPhase` is an integral of `realGradient`, and `realAmplitude`
uses its exponential; neither is asserted rational here.  No global parameter
certificate, compatible `AxisSpace` norm, pressure enclosure, amplitude/phase
enclosure, or paper-exact velocity follows from this input bridge.
