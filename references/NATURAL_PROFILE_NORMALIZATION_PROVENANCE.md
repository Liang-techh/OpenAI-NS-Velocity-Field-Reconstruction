# Natural profile swirl normalization provenance

Status: **formal-structure only; focused implementation validation pending**.

This increment records the pinned rescaling convention and the corresponding
local `NaturalProfileAssembly` correction.  The scaled angular field is the
smooth ratio (F=E/\sqrt{2X}).  It must be passed through the local
`LeadingProfile.F` callback, while the physical swirl field is

```text
F(X, eta) = a(eta) * phi(Lambda * X, eta)
E(X, eta) = sqrt(2 * X) * F(X, eta).
```

The previous direct mapping `LeadingProfile.E = a * phi` was incorrectly
normalized: it gave a nonzero value at the axis and made the velocity bridge
interpret a ratio as physical swirl.  The corrected local assembly exposes
`F`, constructs `E` with `sqrt(2 * X)`, and keeps `paper_exact` false.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- `NavierStokes/NaturalProfile.lean:21-33` defines
  `rescalePoint`, `pullback`, and
  `angularProfile a Lambda Phi = a(eta) * Phi(Lambda * X, eta)`.
- `NavierStokes/NaturalProfile.lean:559-599` defines the final fields
  `ProfileFamily.f`, `U`, `Ubar`, and `Pi`.
- `NavierStokes/ProfileHistories.lean:301-328` defines
  `Profiles.E = sqrt(2 * X) * Profiles.f` and
  `Profiles.pressure = pressure0 + primitive(f^2)`.
- `NavierStokes/NaturalAxisBridge.lean:321-342` gives the scaled axis,
  average, and pressure identities.
- `NavierStokes/NaturalAxisData.lean:27-38` defines
  `D`, `A`, `d`, `L`, `U`, `H`, `W`, and `Z`.
- `NavierStokes/NaturalAxisCoefficients.lean:26-33` fixes the parameter
  window to `(-11/10, 11/10)` and places `[-1, 1]` in its interior.

## Rescaled fields and axis data

With `Y = Lambda * X`, the pinned final fields are

```text
F(X, eta)       = a(eta) * Phi(Y, eta)
E(X, eta)       = sqrt(2 * X) * F(X, eta)
U(X, eta)       = UStar(eta) + Lambda^-1 * u(Y, eta)
radialAverage   = UStar(eta) + Lambda^-1 * B(Y, eta)
Pi(X, eta)      = P0(eta) + Lambda^-1 * P(Y, eta).
```

The scaled fields satisfy `Phi(0, eta) = 1`, `u(0, eta) = 0`,
`B(0, eta) = 0`, and `P(0, eta) = 0`.  Hence the unscaled axis data are
`F(0, eta) = a(eta)`, `E(0, eta) = 0`,
`U(0, eta) = radialAverage(0, eta) = UStar(eta)`, and
`Pi(0, eta) = P0(eta)`.  The source uses `UStar(eta) = 4 * eta + j`.

The unscaled domain is the pullback of
`AxisEvaluation.strip window 20`: `-20 < Lambda * X < 20` and
`-11/10 < eta < 11/10`.  Physical parameter checks commonly use the closed
subinterval `-1 <= eta <= 1`.

## Local bridge requirement

`src/openai_ns_reconstruction/natural_axis.py` must retain the split API:

```text
NaturalProfileAssembly.F  -> smooth ratio a * phi(Lambda * X, eta)
NaturalProfileAssembly.E  -> sqrt(2 * X) * F
to_leading_profile()      -> LeadingProfile(..., F=self.F, E=self.E, ...)
```

This agrees with `velocity.py`, which evaluates physical swirl as
`sqrt(2 * X) * profile.F` when `F` is supplied, and with
`profiles.py`, whose pressure derivative is `F^2`.  Supplying the ratio as
`LeadingProfile.E` without `F` would divide it by `sqrt(2 * X)` and fail at
the regular axis.

## Final pressure boundary

The final profile pressure is the primitive of the squared smooth ratio:

```text
P(Y, eta) = integral_0^Y (a(eta) * Phi(y, eta))^2 dy
Pi(X, eta) - P0(eta) = integral_0^X F(s, eta)^2 ds.
```

This is the `Profiles.pressure` / `ProfileFamily.Pi` path.  It is distinct
from the axial remainder forcing

```text
d * partialEta(P) - 4 * A * eta * P - 2 * eta * Y * partialY(P),
```

which appears only inside `NaturalAxisBridge.axialRemainder`.  A final wide
pressure bridge from the solved `phi` family remains downstream; the landed
coefficient pressure forcing must not be promoted to final `Pi`.

## Truth boundary

Root's focused implementation check was run as
`python -m pytest -q tests/test_natural_axis.py -W error`; it exited `0` with
`4 passed in 0.15 s`.  This documentation task did not rerun that check.

This document records a normalization correction only.  It does not provide
a weighted `AxisSpace` certificate, a converged coefficient-space fixed point,
the final solved `phi/u` family, a final pressure evaluation, or a paper-exact
velocity profile.  The focused result above does not certify any of those
downstream properties.
