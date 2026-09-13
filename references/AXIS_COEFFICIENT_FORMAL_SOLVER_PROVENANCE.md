# Generic triangular formal coefficient solver provenance

Status: **formal-structure only**.

The artifact is `src/openai_ns_reconstruction/axis_coefficient_formal_solver.py`.
Its public binding is `formal_axis_coefficient_solver(x1)`, where `x1` must be
the genuine `ActualScheduleWideFirstPicardState`.  The state exposes
`Lambda`, `epsilon`, `amplitude_log(eta)`, and
`jet_pair(n, m, eta) -> (phi, u)`.  Each returned component is the sparse
`MixedScaleCoefficient` map from `(q,p)` to the normalized numerator of the
channel `a(eta)^q Lambda^-p`.

## Pinned equations and executable binding

The source pin is `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/AxisContraction.lean:327-365` gives the natural remainder,
  including `lin1`, `quad1`, `slow1`, `lin2`, `slow2`, and pressure.
- `NavierStokes/AxisOperators.lean:86-95` gives radial convolution and the
  eta Leibniz product.
- `NavierStokes/AxisOperators.lean:501-504` and `512-518` give the regular
  inverse and differential-family actions.
- `NavierStokes/AxisOperators.lean:681-717` binds the primitive,
  parameter-primitive, Euler, average, multiplication-by-`Y`, and regular
  inverse operators.
- `NavierStokes/AxisWeightEstimates.lean:388` pins the radial divisor
  `n * (n + r - 1)`.
- `NavierStokes/NaturalAxisCoefficients.lean:27` pins the eta window
  `[-11/10, 11/10]`.

The implementation rebuilds fixed fields from
`actual_schedule_axis_coefficient_data(x1.reference)` and takes the amplitude
Bell family from the x1 wide natural source.  It applies

```text
R_phi = naturalResolvent[inverseL * (lin1 + quad1 - slow1/Lambda)]
R_u   = inverseL * (lin2 - slow2/Lambda + pressure)
x     = reference + shift_lambda(R, 1) / 2.
```

The pressure source is the full-derivative normalized family
`a^2 * phi^2`, with
`J1[-4*A*eta*primitive(S) + d*etaDerivative(primitive(S))
    - 2*eta*mulY(S)]`.
The q2 amplitude family uses the existing amplitude-square Bell recurrence;
it is multiplied as a full eta-jet and is never manufactured by a bare
amplitude-key shift.  No caller-fitted coefficient, pressure, epsilon, or
Lambda is accepted.

## Radial causality and stabilization

For a fixed `(eta, n, m)`, all memoization is local to that `jet_pair` call.
The raw operators have the following dependency structure:

1. `J_r` and the primitive read only row `n-1`; `dot`, `param`, and `mixed`
   are products followed by the same radial inverse.
2. The fixed AxisData fields, including `inverseL`, have radial degree zero.
   The pressure chain has the pinned positive radial gain, and its row-zero
   output is structural zero.
3. The angular natural resolvent evaluates the raw row and subtracts only
   previously resolved angular rows, with
   `-4 * reference.phi.jet(1, m, eta)` as its chi family.

Consequently `R(x)[n,m]` reads only `x[i,*]` with `i < n`.  The row-zero
reference is the induction base; the same induction shows that a finite
Picard iterate agrees with the formal coefficient recursion through every row
that has stabilized.  This is a coefficientwise formal stabilization fact,
not a weighted-space convergence theorem.

## Truth boundary

The solver retains `formal_coefficients_materialized = true` and
`fixed_point_materialized = false`, `fixed_point_convergence_certified =
false`, `global_axis_norm_certified = false`, and `paper_exact = false`.
The measured focused check was
`python -m pytest -q tests\\test_axis_coefficient_formal_solver.py -W error`
-> **3 passed in 1.02 s**.  Those tests cover an eta-binomial hand polynomial
and an `a^4` channel, row-zero reference/false-flag behavior, and `n=1,2`
comparisons with the landed `x2`; they do not establish `m=1` or higher-row
coverage.

The solver does not establish the global AxisSpace supremum, weighted-space
membership, contraction convergence, exact paper profile, final velocity,
forcing, or a real-arithmetic proof.

The repository base observed for this addendum is
`67ab5ecd32f80225df8bce2760b943c20a7d3d7b`.  The official Lean source pin is
recorded above; no Lean build or paper-exact promotion is claimed here.
