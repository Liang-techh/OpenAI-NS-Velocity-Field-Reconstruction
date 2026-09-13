# Formal Picard family and radial-filtration provenance

Status: **formal-structure only; theorem-side identification is conditional and
the exact runtime jet bridge remains open**.

The landed `FormalAxisCoefficientSolverState` is a rowwise oracle for the
pinned coefficient formulas.  It is not itself a compatible Banach-space
fixed point.  The exact statement below explains why its row `n` agrees with
the pinned Picard family after at most `n` updates, provided the theorem-side
fixed point and all required analytic hypotheses exist.

## Pinned source

Official source: `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

- `NavierStokes/AxisContraction.lean:339-365` defines the complete
  `naturalRemainder`, including `j1`, `j2`, products, averages, dot/parameter/
  mixed terms, and the pressure branch.
- `NavierStokes/AxisContraction.lean:367-369` defines the reference pair
  `x0 = (S one, -(1/2) * j1(inverseL * zStar))`.
- `NavierStokes/AxisOperators.lean:90-95` gives the product coefficient as a
  finite radial antidiagonal convolution.
- `NavierStokes/AxisOperators.lean:182-190` gives the coefficient projection
  of every bounded linear jet lift; `:380-390` gives the primitive, inverse,
  and `mulY` predecessor-row scales.
- `NavierStokes/AxisOperators.lean:494-508` proves the inverse row-zero and
  predecessor-row identities.
- `NavierStokes/AxisResolvent.lean:321-345` defines the angular filtration and
  proves that the natural linear operator raises radial order.
- `NavierStokes/AxisResolvent.lean:398-405` defines the pinned angular
  operator (Q=(1/2)J_2M_\chi); `:464-471` gives the resolvent equation and
  uniqueness.
- `NavierStokes/AxisReference.lean:65-105` derives the reference coefficient
  recurrence and its factorial closed form from the resolvent equation.
- `NavierStokes/AxisContraction.lean:536-545` gives existence and uniqueness
  of the actual fixed point under the complete-space, norm, amplitude, and
  scale-threshold assumptions.
- `NavierStokes/AxisCoefficientSpace.lean:55-69` gives continuity of coefficient
  evaluations, and `:202-220` gives the weighted evaluation and difference
  bounds.
- `NavierStokes/AxisCoefficientSpace.lean:441-455` gives the derivative bound
  and norm characterization used to turn an all-index coefficient bound into
  an `AxisSpace` witness.
- `NavierStokes/AxisResolvent.lean:312-319` gives the all-index jet bound
  needed to turn a compatible coefficient provider into an `AxisSpace` norm
  witness.

## Radial causality

Write (R(x)) for `naturalRemainder` and

```text
T(x) = x0 + (1 / (2 * Lambda)) * R(x).
```

For a pair of inputs (x,y), suppose every angular and axial coefficient row
below `k` agrees, for every eta-derivative order.  Then every row through `k`
of `R(x)` and `R(y)` agrees.

The reason is structural:

1. A product at row `q` is a finite sum over `i + j = q`, so an outer `j1` or
   `j2` reads only product row `n - 1` when its output row is `n`.
2. `average` preserves radial degree, while `dot1`, `dot2`, `param1`,
   `param2`, `mixed1`, and `mixed2` are inverse bilinear maps with the same
   predecessor-row radial lift.  Thus all `lin`, `quad`, and `slow` branches at
   output row `n` read input rows strictly below `n`.
3. The pressure source is
   `product(product(a,a), product(phi,phi))`.  Its `primitive`,
   `parameterPrimitive`, and `mulY` branches already read predecessor rows,
   and the outer `j1` adds another radial gain.  It therefore also reads only
   rows below `n` (in fact, the source pressure path has gain at least two).
4. The angular resolvent solves `S + Q S = A`.  The pinned `Q` raises radial
   order, so row `n` of `S` is row `n` of `A` minus a finite eta-Leibniz sum of
   resolved rows below `n`.  It cannot reintroduce an input row `n`.
5. Every `j1`/`j2` output row zero is zero.  Hence `R(x)[0,m] = 0` for every
   derivative order and every admissible eta.

The fixed `AxisData` fields are the compatible radial-degree-zero data in
`NaturalAxisBridge.lean:257-269`.  The statement is an exact coefficient
identity assuming the pinned compatible operator definitions; the local
Decimal96 evaluator is only a rounded numerical realization of these maps.

## Picard stabilization

Let (x^{(0)}=x_0) and (x^{(r+1)}=T(x^{(r)})).  Since `R` has zero row zero,
the fixed-point equation gives (x^*[0]=x_0[0]).  Inductively, if

```text
x^(r)[n,m] = x*[n,m]  for every n <= r and every m,
```

then radial causality gives

```text
R(x^(r))[n,m] = R(x*)[n,m]  for every n <= r + 1,
```

and therefore (x^{(r+1)}) agrees with (x^*) through row `r + 1`.
Consequently, the formal row `n` is the row-`n` value of `x^(n)` and remains
unchanged in all later Picard iterates.  This is a finite stabilization claim,
not a convergence rate or a substitute for the global contraction theorem.

The local implementation follows the same triangular equations in
`src/openai_ns_reconstruction/axis_coefficient_formal_solver.py`:
`FormalAxisCoefficientSolverState._build_families` constructs the families;
`get_angular_r` uses the lower resolved row recurrence, while `get_axial_r` is
raw and `get_phi`/`get_u` add the fixed reference row.  `jet_pair` and
`jet_prefix` expose the rows.  The companion
`src/openai_ns_reconstruction/axis_coefficient_picard_family.py` exposes the
`formal_axis_picard_family_state` API.  Reference values are read from the
actual schedule reference anchor; finite `x1`/`x2` tables are not copied into
the recursion.

## Conditional Banach-limit identification

The stabilization statement can be attached to the pinned fixed point without
independently proving a global weighted norm for every finite formal prefix.
Assume the theorem-side contraction hypotheses produce a compatible fixed point
`x*` for the same `T`, and let `x^(r)` be the exact Picard sequence from `x0`.
The contraction construction in `AxisContraction.lean:227-284` is on a complete
space and invokes `ContractingWith.exists_fixedPoint'`; together with the
fixed-point assumptions recorded at `:536-545`, this supplies the norm-limit
route for `x^(r)` (the finite-filtration corollary itself is manual here).

For each fixed `(n,m,eta)`, `AxisCoefficientSpace.abs_jet_sub_le` at
`:216-220` gives, componentwise,

```text
|J[n,m,eta](x^(r)) - J[n,m,eta](x*)|
    <= |weight(epsilon,n,m)| * ||x^(r) - x*||.
```

The right side tends to zero.  For `r >= n`, radial stabilization makes the
left projection equal to the formal row `n` exactly, so its limit is the
coefficient of `x*`.  Thus the exact-arithmetic formal solver row is the
coefficient projection of the compatible fixed point.  The same argument uses
`continuous_jet` at `:55-69`; `abs_jet_le` at `:202-208` bounds the limiting
projection itself.  This is a conditional mathematical identification and is
not currently a separately proved Lean theorem or a Decimal96 equality proof.

## Required `AxisSpace` bridge

The theorem-side route needs one shared `epsilon`, `Lambda`, `AxisData`,
reference pair, and natural remainder.  It must establish:

1. The pinned operator, natural-resolvent, normalized-amplitude, and scale
   hypotheses used by `exists_unique_natural_fixedPoint`.
2. A genuine compatible fixed point `x*` in the complete `AxisSpace` product,
   with its all-index norm witness.  This witness supplies the global bound for
   the formal family through the conditional identification above; a separate
   bound for each finite prefix is unnecessary.
3. The adjacent-eta derivative identity and coefficient projection interface
   required by `AxisCoefficientSpace`, so theorem-side evaluations are the same
   objects as the local solver's `(n,m,eta)` values.
4. An exact or interval-certified runtime bridge from the Decimal/float schedule
   jets to those theorem-side projections.  A finite Decimal comparison is only
   diagnostic until this bridge is supplied.

The current local reference-state adapter intentionally exposes no compatible
fixed-point object or exact jet oracle.  Its formal rows therefore remain a
conditional representation of the theorem-side coefficients.

## Truth and numerical boundary

The filtration and stabilization statement is conditional exact arithmetic.
The local solver uses floating schedule jets and 96-digit Decimal operations;
those operations are not interval arithmetic and do not prove the exact
runtime-to-theorem coefficient equalities, adjacent-eta compatibility, or
paper-exact velocity.  The pinned Lean sources provide the contraction and
evaluation ingredients, but the finite-filtration identification above is a
manual corollary rather than a Lean build result.  No tests were run for this
provenance task.
