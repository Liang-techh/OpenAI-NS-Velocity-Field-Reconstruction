# Conditional formal radial-tail bound provenance

Status: **formal-structure only; conditional analytic estimate**.

Artifacts:

- `src/openai_ns_reconstruction/axis_coefficient_radial_tail.py`
- `references/provenance_manifest_addendum_axis_coefficient_radial_tail.json`

The public function
`axis_coefficient_radial_tail_bound(norm_upper, epsilon, Y, max_n,
eta_order=0, radial_order=0, average=False)` evaluates a conditional majorant
for rows omitted after a finite radial prefix.  `norm_upper` is an explicit
caller-supplied global `AxisCoefficientSpace` norm upper bound.  The function
does not sample coefficient rows, infer a norm from finite data, or attach a
guessed bound to the formal solver.

## Weight input and derived estimate

The pinned formal source is `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/AxisWeightEstimates.lean:24-27` gives

  ```text
  w_epsilon(n,m) = 20^-n epsilon^-m m! binom(n+m,m)
                    / ((n+1)^2 (m+1)^2);
  ```

- `NavierStokes/AxisCoefficientSpace.lean:400-414` supplies the implication
  that a global coefficient-space norm bound `M` controls each coefficient by
  `M * w_epsilon(n,m)`.

Those lines provide the weight and conditional coefficient inequality.  The
radial-tail sum in this module is a derived analytic estimate, not a pinned
Lean tail theorem.  For `N = max_n`, `m = eta_order`, `r = radial_order`, set

```text
L = max(N + 1, r),   J = L - r,   K = m + r,   q = abs(Y)/20.
```

The omitted derivative rows start at `n = L`.  The falling radial factor and
the weight factorial combine as

```text
n!/(n-r)! * m! * binom(n+m,m)
  = K! * binom(n-r+K,K).
```

Bounding the denominator `(n+1)^(2+a)` by `(L+1)^(2+a)`, where
`a = int(average)`, leaves the negative-binomial tail

```text
H = (1-q)^(-(K+1))                                      if J = 0,
H = q^J (1-q)^(-(K+1))
    * sum(k=0..K) binom(J-1+k,k) (1-q)^k                 if J > 0.
```

The returned rational majorant is therefore

```text
M * K! / (20^r epsilon^m (m+1)^2 (L+1)^(2+a)) * H.
```

The implementation handles `q = 0` and `J > 0` as a zero omitted tail, while
`q = 0`, `J = 0` retains the first derivative row.  It requires
`abs(Y) < 20`, finite `M >= 0`, and positive finite `epsilon`.  All supplied
scalar arithmetic is converted to exact `Fraction` values; integer factorials
and binomials remain exact.  Only the final nonnegative rational is converted
to Decimal with 96-digit `ROUND_CEILING`.  Evaluated coefficient roundoff is
not included.

The function has no `eta` argument.  A caller attaching this bound to an
actual solver must separately enforce the solver's pinned eta window
`[-11/10, 11/10]` from `NavierStokes/NaturalAxisCoefficients.lean:27` and
must supply an independently established global norm bound.

## Truth boundary

`AxisCoefficientRadialTailBound` records
`conditional_on_global_norm=True`, `includes_coefficient_roundoff=False`, and
`paper_exact=False`.  These flags are fixed by the result type.  The estimate
does not certify global `AxisSpace` membership, does not prove that the formal
solver's rows satisfy the global norm hypothesis, and does not establish
convergence, a fixed point, a paper-exact velocity, or a physical-coordinate
profile.  The result is a conditional numerical majorant, not exact real or
interval arithmetic after its final Decimal conversion.

The repository base for this addendum is main commit
`67ab5ecd32f80225df8bce2760b943c20a7d3d7b`.
