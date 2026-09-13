# First-Picard `slow2` average/mixed branch provenance

Status: **formal-structure only**.

This increment materializes exactly one remaining mixed-scale branch in the pinned `AxisContraction.naturalRemainder.slow2` expression at the genuine theorem-selected first Picard state `x1`:

```text
product d (mixed1 (average u1) u1)
```

The branch is evaluated on the actual SchedulePressure coefficient data. `d` is the fixed `AxisData.d = 1-eta^2` field from the landed coefficient family, `average` is the pinned radial averaging operator, and `mixed1` is `inverseMixed 1` from `AxisOperators`.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main files:
  - `NavierStokes/AxisContraction.lean`
  - `NavierStokes/AxisOperators.lean`
  - `NavierStokes/NaturalAxisCoefficients.lean`
- Relevant symbols:
  - `NavierStokes.AxisContraction.naturalRemainder`
  - `NavierStokes.AxisOperators.inverseMixed`
  - `NavierStokes.AxisOperators.average`
  - `NavierStokes.NaturalAxisCoefficients.CoefficientFamily.axisData`

The pinned inverse-mixed coefficient identity used here is, for `n >= 1`,

```text
sum_{i+j=n-1} sum_{k+l=m}
  j * binom(m,k) * A[i,k+1] * B[j,l] / n^2.
```

The first argument is `average(u1)` and therefore receives the additional eta derivative. The actual `AxisData.d` field is radial-degree zero, so the outer coefficient product has only the ordinary eta-Leibniz convolution at fixed output radial row.

## Theorem-scale representation

The already-landed first Picard axial state is retained as

```text
u1 = U0 + U1/Lambda + U2/Lambda^2 + a^2 P/Lambda,
```

where `a = exp(Lambda*realPhase)/C` and `Lambda/C` come only from the actual SchedulePressure wide scale chain. The new branch keeps separate:

- ordinary `Lambda^0 .. Lambda^-4` numerators;
- pressure-linear `a^2 Lambda^-1 .. Lambda^-3` numerators;
- pressure-square `a^4 Lambda^-2` numerator.

The `a^2` and `a^4` pieces remain signed-log encoded. No caller-supplied `d`, `Lambda`, `C`, amplitude, pressure state, coefficient table, derivative table, radial cutoff, or surrogate `x1` is accepted.

## Verification

Regression reconstructs every scale numerator independently from the literal inverse-mixed formula followed by the literal eta-Leibniz multiplication by the actual `AxisData.d` jets. The ordinary `Lambda^0` piece is also cross-checked against the previously landed binary64 operator composition

```text
product(d, mixed1(average(reference.u), reference.u)).
```

Additional checks cover exact row-zero vanishing, signed-log `a^2/a^4` scales, failure of binary64 projection for current nonzero theorem-scale pressure pieces, theorem-window/index guards, and rejection of surrogate state objects.

## Deliberate non-claims

This is one `slow2(x1)` branch only. It does **not** materialize the remaining `-param1(u1, d*u1)` branch, complete `slow2(x1)`, complete `naturalRemainder(x1)`, `x2`, a converged/final coefficient-space fixed point, derived final average/pressure fields, `NaturalProfileAssembly`, global all-index weighted `AxisSpace` certification, or the final support/moment/matching/cone checks.

Accordingly Stage 1 remains `formal-structure`, `paper_exact_velocity_available=false`, and `full_reconstruction=false`.
