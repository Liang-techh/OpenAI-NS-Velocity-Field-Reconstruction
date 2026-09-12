# Wide/log naturalRemainder source at the actual reference pair

Status: **formal-structure only**. This artifact does not make the Stage-1 leading profile paper-exact and does not change `paper_exact_velocity_available=false`.

## Pinned source

This increment is pinned to `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, especially `NavierStokes/AxisContraction.lean` and `NavierStokes/NaturalAxisCoefficients.lean`.

The official `AxisContraction.naturalRemainder` defines the pressure source as

`source = O.product (O.product a a) (O.product phi phi)`,

with `t = Lambda^{-1}` and `a = realAmplitude = exp(Lambda*realPhase)/C`. The same file takes the contraction center to be `referencePair O d S`. Therefore, at `x0 = referencePair`, the source is exactly `a^2 * phi0^2`.

## Executable representation

`axis_coefficient_wide_natural_source.py` materializes this one source term without narrowing the theorem-selected amplitude into binary64.

The factory accepts only actual `TailData` and `j`. It rebuilds the landed `ActualScheduleAmplitudeLogState`, inherits the actual SchedulePressure-selected `sigma/epsilon/Lambda/C`, obtains `phi0` from the actual reference pair, and computes `phi0^2` using the already-landed pinned coefficient product. No caller-supplied amplitude, pressure datum, scale, coefficient table, radial cutoff, or derivative table is accepted.

Because `a` has radial degree zero, the radial coefficient row is inherited directly from `phi0^2`. For eta derivative order `m`, the implementation factors out the common nonzero amplitude scale and evaluates

`d_eta^m(a^2 p_n) = a^2 * sum_{k=0}^m binom(m,k) B_k(2 f',...,2 f^(k)) * d_eta^(m-k) p_n`,

where `f = log a`, `f' = Lambda*realGradient`, `p = phi0^2`, and `B_k` is the complete exponential Bell recurrence. The common `a^2` magnitude is stored as `2*log(a)` and the remaining signed Decimal factor is stored separately before conversion to `SignedLogCoefficientJet`.

A legacy `binary64_state()` projection exists only as a fail-closed boundary check. On the current theorem-selected conservative scale, a nonzero source coefficient underflows binary64 and raises `ArithmeticError` rather than becoming zero.

## Validation

Regression coverage checks:

- the source uses the same actual theorem-selected epsilon/reference pair as the amplitude state;
- zeroth eta jets factor exactly as `a^2 * (phi0^2)_n` with the scale kept separate;
- first eta jets satisfy the direct product rule `(a^2 p)'/a^2 = 2 Lambda realGradient p + p'`;
- second eta jets satisfy the explicit Bell/Leibniz identity using `q1=2 Lambda g` and `q2=2 Lambda g'`;
- the current nonzero source refuses binary64 underflow;
- invalid coefficient indices and eta values fail closed.

## Remaining boundary

This is the first wide/log arithmetic bridge inside the genuine `naturalRemainder(x0)` path, but it is not the full remainder. The source still has to be propagated through the pinned linear pressure chain `primitive`, `parameterPrimitive`, `mulY`, coefficient multiplications, and `j1` while retaining its separated amplitude scale. That wide pressure component must then be combined with the ordinary axial remainder without erasing it, the theorem-selected `t=1/Lambda` and fixed-point scale `1/(2 Lambda)` must be represented without binary64 range loss, and only then can the first genuine `naturalRemainder(x0)` / Picard `x1` be evaluated. Global all-index weighted `AxisSpace` certification, fixed-point closure, derived average/pressure profiles, `NaturalProfileAssembly`, and support/moment/matching/cone verification remain open.
