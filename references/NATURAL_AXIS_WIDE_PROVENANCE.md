# Wide natural-axis profile provenance

This addendum records the bounded physical assembly in
`src/openai_ns_reconstruction/natural_axis_wide.py`. The implementation
consumes the formal sparse coefficient graph and returns
`wide_natural_profile_prefix(solver, max_n, X, eta)` as a frozen
`WideNaturalProfilePrefix`.

The solver's radial coordinate is `Y = Lambda * X`. The public physical chart
accepts finite `Decimal X` with `X >= 0`, computes the finite Decimal product
`Y` exactly at widened precision, and requires `0 <= Y < 20` and `|eta| <= 1`.
The maps retain `(q, p)` channels for `a^q Lambda^(-p)`. Odd amplitude powers
are allowed because the physical angular ratio is `F = a phi`; no amplitude
power is formed and no sparse map is collapsed to binary64.

The assembled maps are:

* `F`: the value angular jet multiplied by the symbolic `q = 1` channel;
* `E`: `sqrt(2 X) F`, with exact zero on the axis;
* `U`: the ordinary `4 eta + j` term plus the axial jet shifted by one
  inverse-Lambda channel;
* `dU_deta`: ordinary `4` plus the `eta_order = 1` axial jet shifted by one
  inverse-Lambda channel;
* `average_U`: ordinary `4 eta + j` plus the shifted radial average;
* `d_average_U_deta`: ordinary `4` plus the shifted eta derivative of the
  radial average;
* `V0`: the sparse regular radial-velocity numerator
  `X/L * (2 eta U - 2 D eta average_U - (1 - eta^2) d_average_U_deta)`, with
  structural zero at `X = 0`;
* `Pi`: the schedule axis-pressure quadrature value plus the primitive pressure
  prefix shifted by one inverse-Lambda channel.

The transport geometry uses the actual `solver.data.D` and `solver.data.h`
values imported with `Decimal.from_float`.  Evaluation validates `|eta| <= 1`
and `L = 1 - 2 h eta^2 > 0`; the sparse products and scalar factors use a
96-digit local Decimal context.

`terms_log(field_name)` returns signed-log terms and preserves
`q * amplitude_log` separately from `ln(abs(numerator)) - p ln(Lambda)`.
`conditional_truncation_bounds(budget)` returns an immutable eight-field mapping
whose values are upward-rounded Decimal omitted-row bounds. The new
`d_average_U_deta` tail delegates to the eta-order-one averaged bound. The
`V0` tail combines the three existing velocity tails with the exact rational
factors from `X`, `eta`, `D`, `h`, and positive `L`, then rounds the final result
upward to Decimal96. It is explicitly
conditional on identifying the supplied actual budget with the global
AxisSpace norm; coefficient roundoff and schedule axis-pressure quadrature
roundoff are excluded.
`paper_exact`, global axis-norm certification, fixed-point materialization,
truncation certification, and coefficient-roundoff certification remain false.
The schedule axis-pressure term is a quadrature approximation; it is not the
pressure forcing used by the remainder graph and contributes no certified
quadrature error to the conditional bounds.

The source formulas are pinned to `NaturalProfile.lean` lines 559--599,
`ProfileHistories.lean` lines 301--328, and `NaturalAxisBridge.lean` lines
321--342 at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`. This artifact does not establish
global `AxisSpace` membership, convergence, interval error budgets, or a
paper-exact reconstruction.
