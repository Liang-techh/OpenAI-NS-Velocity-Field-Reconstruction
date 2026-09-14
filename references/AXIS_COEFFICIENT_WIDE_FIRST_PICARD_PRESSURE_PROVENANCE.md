# Axis coefficient wide first-Picard pressure provenance

## Scope

This increment materializes only the pinned pressure constituent of
`naturalRemainder(x1)` from the genuine first-Picard source already materialized
by PR #323:

`source(x1) = a^2 * phi1^2`.

It does not rebuild the phase/reference hierarchy and does not accept a
caller-supplied source, pressure table, amplitude, Lambda, C, or coefficient
family.

## Pinned source

Formal source repository: `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

`NavierStokes.AxisContraction.naturalRemainder` defines

`pressure = j1 (-(4*A*eta)*primitive(source)
                + d*parameterPrimitive(source)
                -(2*eta)*mulY(source))`.

`NavierStokes.AxisOperators` implements `primitive`,
`parameterPrimitive`, `mulY`, and `j1 = regularInverse(...,1)` as genuine
compatible coefficient-space operators. Their row formulas are replayed
directly here: `primitive` and `parameterPrimitive` use predecessor row divided
by `n`, `mulY` uses the predecessor row, and `j1` uses predecessor row divided
by `n^2`; row zero is zero by the pinned definitions.

## Executable realization

`axis_coefficient_wide_first_picard_pressure.py` consumes only
`ActualScheduleWideFirstPicardNaturalSourceState`. The source's five
`Lambda^0..Lambda^-4` normalized families are propagated independently through
the linear pressure chain. Multiplication by ordinary `eta` and `d` coefficient
fields uses the complete finite radial convolution and eta-Leibniz sum. Those
ordinary fields come from the same actual-schedule `AxisData` already bound into
x1.

The common nonzero `a^2` magnitude is never formed in binary64. Each pressure
coefficient jet keeps five Decimal numerators plus the original amplitude log;
`pressure_terms_log()` restores each term as
`a^2 * p_k / Lambda^k` in split signed-log coordinates.

The pinned zero output rows are evaluated only after touching the corresponding
upstream source/input coordinate, so invalid upstream data cannot be silently
hidden by a structural zero. No finite radial cutoff or absent coefficient is
interpreted as zero.

## Regression boundary

The regression checks, without tolerance:

- the `primitive`, `parameterPrimitive`, and `mulY` row identities for all five
  inverse-Lambda scales;
- the zeroth-eta pressure-input formula against the independently exposed
  actual `eta`, `d`, and `A` fields;
- the final `j1` predecessor-row / `n^2` identity and its pinned row-zero rule;
- preservation of the common signed-log `a^2` scale and each explicit
  inverse-Lambda denominator;
- fail-closed rejection of surrogate source states and invalid coordinates.

This establishes only `pressure(x1)`. Complete `naturalRemainder(x1)` still
requires the remaining angular `lin1/quad1/slow1` and axial `lin2/slow2`
constituents to be evaluated on x1 and recombined with this pressure. Therefore
x2, Picard convergence, all-index `AxisCoefficientSpace` membership, fixed-point
phi/u, derived average/pressure fields, `NaturalProfileAssembly`, and
paper-exact velocity remain unproved/unmaterialized.
