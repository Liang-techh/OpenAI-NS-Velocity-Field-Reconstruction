# Axis coefficient product provenance

Status: **formal-structure**. This increment does **not** make the leading profile or velocity paper-exact.

## Pinned official source

Official Lean repository: `openai/NavierStokesAndEuler` at commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The executable adapter in `src/openai_ns_reconstruction/axis_coefficient_product.py` follows these pinned definitions:

- `NavierStokes/AxisWeightEstimates.lean`
  - `jetProduct`: for radial degree `n` and parameter derivative order `m`, sum over `i+j=n` and `k+l=m` with the exact binomial coefficient `m.choose k`.
- `NavierStokes/AxisOperators.lean`
  - `leibnizSum`: the finite parameter-derivative Leibniz sum.
  - `productFamily`: radial convolution of those actual Leibniz sums.
  - `product`: the bounded bilinear operator obtained from `productFamily` after compatibility and weighted-norm estimates are proved.
  - `productFamily_bound`: the official all-index estimate uses the constant `64` once genuine `AxisSpace` norms are available.

The Python product is intentionally landed at the `productFamily` level. It evaluates the actual unnormalized jets already carried by `AxisCoefficientJetState`, and it requires the two inputs to have exactly the same theorem-selected `epsilon`. It accepts no replacement coefficient table, no fitted derivative data, and no independent weight scale.

## Executable identity

For operands `A,B`, the returned lazy state uses

```text
J_AB[n,m](eta)
  = sum_{i+j=n} sum_{k+l=m}
      binom(m,k) J_A[i,k](eta) J_B[j,l](eta).
```

The radial and parameter sums are finite for each requested `(n,m)`. Every factor is an actual parameter derivative supplied by the operand state. Non-finite intermediate products fail closed instead of becoming stored coefficient values.

## Independent regression boundary

`tests/test_axis_coefficient_product.py` uses the actual SchedulePressure-derived reference states, not a caller-provided surrogate table. It checks:

1. the axial-reference square against independently written first/second derivative identities for `u_1(eta)^2` and the exact radial support consequence that only degree two survives;
2. a mixed angular/axial first parameter jet against a centered finite-difference derivative of the independently assembled zeroth-jet radial product;
3. zeroth parameter jets against a direct radial polynomial convolution;
4. fail-closed rejection of mismatched coefficient-space weight scales and preservation of the existing index/window guards.

The finite-difference comparison is only numerical cross-validation of the executable derivative identity. It is not used to generate production jets and is not a substitute for the pinned Lean compatibility proof.

## What remains unresolved

This increment is one genuine primitive of the official operator algebra, but it deliberately does **not** claim the complete `AxisOperators` / `coefficientOperators` backend. Remaining Stage-1 blockers include:

1. the global all-index weighted norm certificate needed to instantiate genuine `AxisSpace` membership for the materialized reference states and derived states;
2. the remaining pinned linear operators needed by `coefficientOperators` (in particular the regular radial inverse/shift/derivative/average paths and the bounded natural resolvent connection);
3. the official `naturalRemainder` evaluation at the actual reference pair;
4. the first genuine Picard iterate and converged coefficient-space fixed point `phi/u`;
5. derived average/pressure fields and `NaturalProfileAssembly`;
6. support, moments, matching, and cone-condition closure for Theorem 4.6.

The official `productFamily_bound <= 64 ||A|| ||B|| * weight` is recorded here as provenance but is **not** converted into a Python `AxisSpace` membership claim until the missing global input norms are constructively certified. Accordingly `paper_exact` and `paper_exact_velocity_available` remain fail-closed (`false`).
