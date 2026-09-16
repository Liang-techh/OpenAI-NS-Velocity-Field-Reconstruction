# Finite formal mixed-scale profile-prefix provenance

Status: **formal-structure only**.

Artifacts:

- `src/openai_ns_reconstruction/axis_coefficient_formal_solver.py`
- `src/openai_ns_reconstruction/axis_coefficient_profile_prefix.py`
- `references/provenance_manifest_addendum_axis_coefficient_profile_prefix.json`

The solver extension `jet_prefix(max_n, m, eta)` evaluates rows `0..max_n`
from one eta-local recursive family graph.  The public factory
`formal_axis_profile_prefix(solver, max_n, Y, eta, radial_order=0,
eta_order=0)` then evaluates a finite natural-radial polynomial prefix.  `Y`
is a finite nonnegative `Decimal`; no floating-point theorem scale or physical
coordinate transform is introduced.

For a finite sequence of coefficient maps `c_n`, the angular and axial maps
use

```text
P_r(Y) = sum(n >= r) n!/(n-r)! * Y^(n-r) * c_n.
```

The axial-average map applies the pinned radial average to each original row
before this radial derivative/evaluation:

```text
A_r(Y) = sum(n >= r) n!/(n-r)! * Y^(n-r) * c_n/(n+1).
```

At `Y = 0`, only row `n = r` contributes.  Every sparse `(q,p)` channel is
preserved in the `MixedScaleCoefficient` result.  The implementation forms
neither `a(eta)^q` nor `Lambda^-p`.  Signed-log views keep each channel
separate, with `log_scale = q * amplitude_log` and
`log_factor = ln(abs(numerator)) - p * ln(Lambda)`; `q = 0` is valid.

## Pinned source map

The formal source is `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/AxisOperators.lean:86-95` supplies radial and eta-binomial
  coefficient products;
- `NavierStokes/AxisOperators.lean:501-504` supplies the regular inverse
  family used by the coefficient solver;
- `NavierStokes/AxisOperators.lean:512-518` supplies differential-family
  actions;
- `NavierStokes/AxisOperators.lean:681-717` binds radial average and related
  coefficient operators;
- `NavierStokes/AxisWeightEstimates.lean:388` supplies the radial divisor
  convention used by the inverse families;
- `NavierStokes/NaturalAxisCoefficients.lean:27` pins the eta window
  `[-11/10, 11/10]`.

The local `axis_coefficient_average.py` records the row map
`c_n -> c_n/(n+1)`, while `axis_coefficient_mixed_scale.py` supplies the
immutable sparse channel arithmetic.  The prefix is anchored only to the
formal solver's genuine `ActualScheduleWideFirstPicardState`; it does not
replace that solver's recurrence with a finite iterate.

## Truth and numerical boundary

`FormalAxisProfilePrefix` exposes `paper_exact`,
`global_axis_norm_certified`, `fixed_point_materialized`, and
`truncation_certified`, all false.  A finite prefix is not a converged profile,
and the uncomputed tail has no inferred bound.  The artifact does not claim a
global weighted `AxisSpace` estimate, fixed-point convergence, paper-exact
velocity, forcing, or a physical-coordinate evaluation.

Channel arithmetic and prefix scalars are rounded in a local 96-digit
`Decimal` context.  A signed-log `q * amplitude_log` product widens its local
context to retain the full finite Decimal amplitude log before storing the
separate scale component.  This is finite numerical evaluation, not exact
real arithmetic, interval arithmetic, or a theorem certificate.

The repository base for this addendum is main commit
`67ab5ecd32f80225df8bce2760b943c20a7d3d7b`.  Validation is limited to the
bounded `py_compile` check and one small `Y = 0` prefix smoke; no full test
suite is claimed here.
