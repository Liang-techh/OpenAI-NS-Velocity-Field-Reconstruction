# Stage 1 wide-profile local integration report

Status: **reviewed conditional wide-profile increment**.

This report records the local integration boundary for the Stage 1 profile
path.  The wide sparse profile maps and actual-chain conditional F/P budget are
implemented and reviewed.  They remain formal conditional interfaces rather
than a global compatible `AxisSpace` certificate.

## Implemented pressure prefix

The formal coefficient solver now exposes
`profile_jet_prefix(max_n, m, eta)`, returning one eta-local graph of
`(phi, u, P)` sparse mixed-scale rows.  The pressure family is

```text
P = primitive((a * phi)^2) = primitive(a^2 * phi^2).
```

It uses the actual x1 amplitude source's full Bell derivative family and the
generic zero-axis primitive.  Thus `P[0,m]` is structural zero and
`P[n,m] = S[n-1,m]/n` for `n > 0`.  This path is separate from the axial
natural-remainder pressure forcing `J1[d P_eta - 4 A eta P - 2 eta Y P_Y]`;
it applies no outer inverse-Lambda factor and does not assemble physical
`Pi`.

`formal_axis_profile_prefix(...)` now carries a sparse `pressure` map and a
`pressure_terms_log()` view, applying the same finite radial and eta-order
evaluation used for angular and axial maps.  Existing truth flags remain false.

## Wide profile and conditional budget

`natural_axis_wide.py` exposes the reviewed `wide_natural_profile_prefix`,
retaining the actual Decimal Lambda and symbolic amplitude in sparse maps
`F`, `E`, `U`, `d_eta U`, `Ubar`, and `Pi`.  The physical `Pi` assembly is now
present as a formal sparse prefix; it does not claim a converged or paper-exact
profile.

The reviewed actual-chain conditional F/P budget supplies the six-field
`conditional_truncation_bounds` mapping.  It uses exact Fraction scaling and an
upward square-root enclosure, while retaining explicit exclusions for global
AxisSpace identification, coefficient roundoff, and axis-pressure quadrature
roundoff.

## Local evidence

The bounded source check passed:

```text
python -B -m py_compile \
  src/openai_ns_reconstruction/axis_coefficient_formal_solver.py \
  src/openai_ns_reconstruction/axis_coefficient_profile_prefix.py
```

The one actual-schedule pressure smoke passed.  For the regression schedule,
`profile_jet_prefix(1, 0, 0.17)` returned an empty row-zero pressure support and
row one support `((2, 0),)` with numerator `1`, i.e. `P_Y(0)=a^2` at the axis.

The focused profile acceptance command was:

```text
python -m pytest -q tests\test_natural_axis_wide.py tests\test_axis_coefficient_profile_prefix.py -W error
```

It returned `5 passed in 0.66 s`.  The final wide-only acceptance command was:

```text
python -m pytest -q tests\test_natural_axis_wide.py -W error
```

It returned `4 passed in 0.64 s`, covering exact Fraction scalings, upward
square-root enclosure, axis-zero tails, and solver mismatch rejection.  No
full-suite rerun was performed for this documentation increment.

## Remaining boundary

Global weighted `AxisSpace` membership, adjacent-eta compatibility, coefficient
and phase-quadrature roundoff budgets, fixed-point convergence, and paper-exact
certification remain open.  The scaled formal pressure primitive and sparse
physical `Pi` assembly are implemented conditional prefixes; they are not a
completed fixed point or paper-exact profile.

Repository base: `67ab5ecd32f80225df8bce2760b943c20a7d3d7b`.
