# Wide/log naturalRemainder pressure at the actual reference pair

Status: **formal-structure only**. This artifact does not make the Stage-1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

This increment is pinned to `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, especially `NavierStokes/AxisContraction.lean`, `NavierStokes/AxisOperators.lean`, and `NavierStokes/AxisWeightEstimates.lean`.

The official `AxisContraction.naturalRemainder` defines

`source = a^2 * phi^2`

and then

`pressure = j1 (-(4 A eta) * primitive(source) + d * parameterPrimitive(source) - (2 eta) * mulY(source))`.

At the contraction center `x0 = referencePair`, the landed wide source is the genuine theorem-selected `a^2 * phi0^2` built from the actual SchedulePressure `Lambda/C` chain.

## Executable representation

`axis_coefficient_wide_natural_pressure.py` propagates that exact wide source through the complete pinned linear pressure chain while retaining the same common nonzero `a(eta)^2` scale.

For each coefficient jet, the implementation stores only the Decimal factor `F[n,m](eta)` defined by

`d_eta^m pressure_n(eta) = a(eta)^2 * F[n,m](eta)`.

The individual pinned maps are implemented coefficientwise:

- `primitive`: row predecessor divided by `n`;
- `parameterPrimitive`: row predecessor at eta-jet order `m+1`, divided by `n`;
- `mulY`: row predecessor with no radial divisor;
- multiplication by ordinary `eta` and `d` fields: exact finite radial convolution plus eta Leibniz sum;
- final `j1`: predecessor row divided by `radialDivisor(1,n-1)=n^2`.

No caller-supplied amplitude, pressure datum, `sigma`, `epsilon`, `Lambda`, `C`, derivative table, coefficient table, radial cutoff, or scale factor is accepted. The factory accepts only actual `TailData` and `j`, rebuilds the landed theorem chain, and inherits its actual reference state and coefficient scale.

`jet_log` restores the common `2*log(a)` scale only at the final signed-log boundary. `binary64_state()` is explicitly fail-closed: a nonzero coefficient may either be representable as a finite nonzero binary64 value or raise `ArithmeticError`; it is never silently changed to zero or infinity.

## Validation

Regression coverage checks:

- exact binding to the landed actual wide natural source and coefficient scale;
- exact radial row identities for `primitive`, `parameterPrimitive`, and `mulY`;
- the complete pressure input against an independently written literal Decimal implementation of the pinned formula;
- the final `j1` predecessor-row map and `n^2` divisor;
- nonzero signed-log pressure coefficients cannot silently project to binary64 zero;
- invalid indices and eta-window violations fail closed.

## Remaining boundary

The theorem-selected wide pressure contribution is now materialized at `x0`, but this is not yet the complete `naturalRemainder(x0)`. The next mixed-scale seam is to combine this `a^2`-scaled pressure with the ordinary axial `lin2 - t*slow2` contribution without erasing the wide term, while also representing theorem-selected `t=1/Lambda` and the outer Picard scale `1/(2 Lambda)` without binary64 loss. Only after that combination can the first genuine Picard `x1` be materialized. Global all-index weighted `AxisSpace` certification, fixed-point closure, derived average/pressure profiles, `NaturalProfileAssembly`, and support/moment/matching/cone verification remain open.
