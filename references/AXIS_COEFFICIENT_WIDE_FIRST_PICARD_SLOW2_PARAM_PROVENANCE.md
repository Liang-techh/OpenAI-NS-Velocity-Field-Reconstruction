# First-Picard `slow2` parameter branch provenance

Status: **formal-structure only**.

This increment materializes exactly the remaining parameter-product branch in the pinned `AxisContraction.naturalRemainder.slow2` expression at the genuine theorem-selected first Picard state `x1`:

```text
- param1 u1 (product d u1)
```

The branch is evaluated on the actual SchedulePressure coefficient data. `d` is the fixed `AxisData.d = 1-eta^2` field, and `param1` is `inverseParamProduct 1` from the pinned `AxisOperators` record. No replacement coefficient data are accepted.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Main files:
  - `NavierStokes/AxisContraction.lean`
  - `NavierStokes/AxisOperators.lean`
  - `NavierStokes/NaturalAxisCoefficients.lean`
- Relevant symbols:
  - `NavierStokes.AxisContraction.naturalRemainder`
  - `NavierStokes.AxisOperators.inverseParamProduct`
  - `NavierStokes.AxisOperators.product`
  - `NavierStokes.NaturalAxisCoefficients.CoefficientFamily.axisData`

For `param1 = inverseParamProduct 1`, the pinned coefficient identity is, for output row `n >= 1`,

```text
sum_{i+j=n-1} sum_{k+l=m}
  binom(m,k) * U[i,k+1] * (d*U)[j,l] / n^2.
```

The implementation includes the literal minus sign appearing in `slow2`. Because `AxisData.d` has radial degree zero, `(d*U)[j,l]` uses only the exact eta-Leibniz convolution.

## Theorem-scale representation

The landed first Picard axial state remains

```text
u1 = U0 + U1/Lambda + U2/Lambda^2 + a^2 P/Lambda,
```

where `a = exp(Lambda*realPhase)/C`, and `Lambda/C` come only from the actual SchedulePressure wide-scale chain. The new branch keeps separate:

- ordinary `Lambda^0 .. Lambda^-4` numerators;
- pressure-linear `a^2 Lambda^-1 .. Lambda^-3` numerators;
- pressure-square `a^4 Lambda^-2` numerator.

The `a^2` and `a^4` pieces remain signed-log encoded. No caller-supplied `d`, `Lambda`, `C`, amplitude, pressure state, coefficient table, derivative table, radial cutoff, or surrogate `x1` is accepted.

## Verification

Regression independently reconstructs the literal `product(d,u1)` eta-Leibniz convolution and then the complete `inverseParamProduct 1` coefficient sum at every retained scale. The ordinary `Lambda^0` piece is separately cross-checked against the already-landed binary64 operator composition

```text
- param1(reference.u, product(AxisData.d, reference.u)).
```

Additional checks cover exact row-zero vanishing, signed-log `a^2/a^4` scales, binary64 fail-closed behavior for nonzero theorem-scale pressure pieces, theorem-window/index guards, and rejection of surrogate state objects.

## Deliberate non-claims

All four constituent `slow2(x1)` branches now have individual materializations, but this increment does **not** yet sum them into one typed complete `slow2(x1)` object. It therefore does not claim complete `slow2(x1)`, complete `naturalRemainder(x1)`, `x2`, a converged/final coefficient-space fixed point, derived final average/pressure fields, `NaturalProfileAssembly`, global all-index weighted `AxisSpace` certification, or the final support/moment/matching/cone checks.

Accordingly Stage 1 remains `formal-structure`, `paper_exact_velocity_available=false`, and `full_reconstruction=false`.
