# Complete first-Picard `slow2(x1)` composition provenance

Status: **formal-structure only**.

This increment closes the typed composition seam left after the four pinned
`AxisContraction.naturalRemainder.slow2` constituents were independently
materialized at the genuine theorem-selected first Picard state `x1`.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Primary source: `NavierStokes/AxisContraction.lean`
- Relevant symbol: `NavierStokes.AxisContraction.naturalRemainder`

The composed expression is exactly

```text
j1 ((2*A*eta) * (u1*u1))
+ dot1 ((2*D*eta) * average(u1)) u1
+ d * mixed1 (average(u1)) u1
- param1 u1 (d*u1).
```

## Actual-data binding

The constructor accepts only one genuine
`ActualScheduleWideFirstPicardState`. Every constituent state is recreated from
that identical object and must preserve its theorem-selected `Lambda` and
coefficient-space `epsilon`. At each requested `(n,m,eta)` all four branch jets
must also share exactly the same `Lambda` and signed-log amplitude value before
addition is permitted.

No caller-supplied coefficient table, pressure state, `Lambda`, `C`, amplitude,
cutoff, surrogate `x1`, or replacement branch value is accepted.

## Scale-preserving composition

Each landed branch already separates the same scale basis:

- ordinary `Lambda^0 ... Lambda^-4` numerators;
- pressure-linear `a^2 Lambda^-1 ... Lambda^-3` numerators;
- pressure-square `a^4 Lambda^-2` numerator.

The complete `slow2` state adds only equal-basis Decimal numerators. It does not
convert signed-log pressure scales to binary64 and therefore does not silently
drop nonzero theorem-scale terms by underflow. The implementation performs the
wide composition at the repository's 96-digit Decimal precision. Regression
reconstructs the expected branch sums and amplitude powers under the same
explicit precision rather than relying on the process-global Decimal default.

## Verification boundary

Regression verifies that the complete jet is the precision-consistent exact
componentwise sum of all four landed branch jets, including the exact row-zero
result, shared `x1/Lambda/epsilon/amplitude` identity, signed-log scale retention,
invalid index/window rejection, and surrogate-state rejection.

This increment establishes **complete `slow2(x1)` only**. It does not establish
complete `naturalRemainder(x1)`: the remaining angular `lin/quad/slow` terms and
the source/pressure/resolvent recombination still require theorem-scale
materialization. Consequently it does not materialize `x2`, a converged fixed
point, final `phi/u`, derived average/pressure, `NaturalProfileAssembly`, global
weighted `AxisSpace` certification, or paper-exact velocity.

`full_reconstruction=false` and `paper_exact_velocity_available=false` remain
mandatory.
