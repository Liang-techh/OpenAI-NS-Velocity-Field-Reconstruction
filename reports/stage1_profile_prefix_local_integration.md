# Stage 1 profile-prefix local integration report

Status: **locally reviewed formal-structure increment; conditional radial-tail bound accepted**.

This report records the finite profile-prefix bridge from the generic mixed-scale
solver to radial polynomial evaluation. It does not claim a converged profile,
global weighted coefficient state, final pressure, or a fixed point.

## Repository identity

- Baseline main commit: `67ab5ecd32f80225df8bce2760b943c20a7d3d7b`
- Local integration branch: `codex/stage1-complete-slow2`
- Current public issue scopes remain #1 Stage 1, #2 Stage 2, #3 Stages 3–6,
  and #4 Stages 7–8.

## Measured evidence

The profile-prefix focused target was run as:

```text
python -m pytest -q tests/test_axis_coefficient_profile_prefix.py -W error
```

It exited `0` with `2 passed` in `0.37 s`.

The directed-rounding scalar target was run as:

```text
python -m pytest -q tests/test_axis_fixed_point_picard.py -W error
```

It exited `0` with `7 passed` in `0.17 s`. The test covers the lower-rounded
`2 * Lambda` denominator and lower-rounded `1 - q` denominator against exact
`Fraction` reciprocals, including the `q = 1e-100` and 96-digit `Lambda`
adversarial cases.

The conditional radial-tail and actual profile-budget target was run as:

```text
python -m pytest -q tests\test_axis_coefficient_radial_tail.py -W error
```

It exited `0` with `4 passed` in `0.17 s`. The target includes the Decimal
precision-invariance regression and an actual-parameter-derived norm budget;
it does not claim global `AxisSpace` membership.

The natural-axis normalization target was run as:

```text
python -m pytest -q tests/test_natural_axis.py -W error
```

It exited `0` with `4 passed` in `0.15 s`. `NaturalProfileAssembly` now
exposes `F = a phi` and computes the physical swirl as `E = sqrt(2 X) F`,
including the regular axis value.

No full-suite rerun, remote publication, or commit is recorded here.

## Integrated scope

`FormalAxisCoefficientSolverState.jet_prefix(max_n, m, eta)` supplies one
eta-local graph of rows. `formal_axis_profile_prefix(solver, max_n, Y, eta,
radial_order, eta_order)` evaluates the finite natural-radial polynomial

```text
sum(n >= r) n!/(n-r)! * Y^(n-r) * c_n,
```

and its axial radial-average variant divides each original row by `n + 1`
before the radial derivative/evaluation. Every sparse `(amplitude_power,
inverse_Lambda_power)` channel remains separate. Signed-log views preserve
amplitude and Lambda scales without forming either wide power.

The profile prefix is a finite `Y`-radial formal evaluation. Its truth flags
remain false for paper exactness, fixed-point materialization, truncation
certification, and global weighted `AxisSpace` certification. The conditional
`actual_schedule_profile_budget(solver)` derives its norm upper bound from the
same actual reference norm and scalar one-step certificate; it does not identify
the formal state with a compatible global fixed point or absorb coefficient
roundoff. The pressure forcing carried inside the axial remainder is not the
final profile pressure.

## Remaining boundary

`axis_coefficient_radial_tail.py` now has a focused accepted conditional
majorant under the pinned `AxisSpace` weight, with
`actual_schedule_profile_budget(solver)` deriving its norm upper bound from the
same reference norm and scalar schedule certificate. This remains conditional:
the compatible global `AxisSpace` identification and coefficient-roundoff budget
are still open, so the finite prefix and tail bound do not establish a global
fixed point.

The actual remaining Stage 1 blockers are compatible weighted-space/global
majorants and coefficient-roundoff control, the final scaled pressure
`P = primitive((a phi)^2)` distinct from remainder forcing, and a wide physical
`NaturalProfileAssembly` capable of retaining the actual `Lambda/C` representation.
No complete reconstruction or fixed point is claimed.
