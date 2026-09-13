# Formal scaled-pressure profile bridge provenance

Status: **formal-structure only**.

Artifacts:

- `src/openai_ns_reconstruction/axis_coefficient_formal_solver.py`
- `src/openai_ns_reconstruction/axis_coefficient_profile_prefix.py`
- `references/provenance_manifest_addendum_axis_coefficient_profile_pressure.json`

The formal solver now exposes
`profile_jet_prefix(max_n, m, eta)`, returning one eta-local cached graph of
`(phi, u, P)` sparse mixed-scale rows.  The pressure family is the final
scaled profile pressure

```text
P = primitive((a * phi)^2) = primitive(a^2 * phi^2).
```

It is kept separate from the axial natural-remainder forcing
`J1[d P_eta - 4 A eta P - 2 eta Y P_Y]`.  No outer `1/Lambda` factor or
physical pressure assembly is introduced here.

## Pinned source and executable binding

The formal source is `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/NaturalAxisBridge.lean:321-342` pins the scaled pressure
  primitive used by the natural profile bridge;
- `NavierStokes/NaturalProfile.lean:559-599` binds the profile pressure
  contribution and its radial primitive semantics.

The executable family constructs a full eta-derivative Bell family for `a^2`
from the actual x1 wide natural source helper, multiplies it with the full
eta-jet product `phi * phi`, and applies the generic zero-axis primitive:

```text
S = product(a^2, product(phi, phi))
P = primitive(S)
P[0,m] = 0,
P[n,m] = S[n-1,m] / n  for n > 0.
```

The amplitude factor is therefore not a bare channel-key shift.  Its Bell
derivatives remain part of the eta Leibniz product, and all nonzero `(q,p)`
channels are retained at Decimal96 precision.  The solver continues to use
the genuine x1 reference anchor and formal triangular recursion; finite x1
coefficient totals are not substituted into the pressure family.

`formal_axis_profile_prefix(...)` now obtains `(phi,u,P)` rows from one graph
and applies the shared finite radial evaluator to all three maps.  The frozen
prefix object carries a sparse `pressure` map and `pressure_terms_log()` view.
The radial order and eta order apply to pressure exactly as they do to angular
and axial maps.

## Truth and numerical boundary

The row-zero pressure is structural zero, while the first radial pressure row
is the actual `(2,0)` `a^2` channel from the x1 anchor when the axis angular
reference coefficient is one.  This is a scaled formal coefficient family;
it does not materialize physical `Pi = P0 + P/Lambda`, final velocity, or a
converged fixed point.

The existing prefix flags remain false for `paper_exact`,
`global_axis_norm_certified`, `fixed_point_materialized`, and
`truncation_certified`.  The module does not prove the unknown radial tail,
global compatible `AxisSpace` membership, or any pressure norm estimate.
Decimal96 channel arithmetic is rounded numerical evaluation, not exact real
arithmetic or a theorem certificate.

Validation for this increment is limited to a bounded `py_compile` check and
one actual-schedule smoke confirming `P[0] = 0` and the axis derivative
`P_Y(0) = a^2` in channel `(2,0)` with numerator one.  No full test suite is
claimed here.

The repository base for this addendum is main commit
`67ab5ecd32f80225df8bce2760b943c20a7d3d7b`.
