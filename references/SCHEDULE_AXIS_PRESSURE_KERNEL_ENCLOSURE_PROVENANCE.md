# Strict SchedulePressure kernel and finite eta-jet enclosures

Status: rigorous real scalar/jet enclosure for the pressure kernel only.
Stage 1, the all-real-line pressure integral, and full reconstruction remain incomplete.

## Pinned formula and source boundary

The runtime formula already used by
`src/openai_ns_reconstruction/schedule_axis_pressure.py` and
`src/openai_ns_reconstruction/schedule_axis_pressure_jets.py` is

```text
k_a(eta) = (1 + eta^2)^(-2 a),     0 <= a <= 1.
```

The source mapping is the pinned OpenAI formal revision
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`,
especially `NavierStokes/PressureDatum.lean` and
`NavierStokes/SchedulePressure.lean`. The existing
`pressure_kernel_normalized_taylor` evaluates arbitrary finite real eta jets
with floating truncated-series algebra; it explicitly does not provide a
rigorous interval enclosure. NS006 adds that strict finite-real entrypoint and
does not alter the existing floating evaluator.

## Exact-rational value enclosure

Inputs are exact integers/Fractions or closed `RationalInterval` objects.
Floats, Decimals, complex numbers, out-of-range `a`, and malformed caps fail
closed. For an eta interval, `1+eta^2` is enclosed by exact rational interval
arithmetic. For a general exponent interval the value is written as

```text
k_a(eta) = exp(-2 a log(1+eta^2)).
```

`log(x)` for rational `x>=1` is enclosed with the positive atanh series

```text
log(x) = 2 sum_{j>=0} z^(2j+1)/(2j+1),  z=(x-1)/(x+1),
```

after exact power-of-two range reduction. The omitted positive tail is bounded
by a geometric denominator bound. `exp(-y)` for rational `y>=0` is range
reduced to `q<=1/2`, bracketed by successive even/odd alternating Taylor sums,
then restored by exact repeated squaring. Term/range-reduction caps are
explicit failures, never silent fallbacks.

Two exact real identities bypass transcendental enclosure entirely:

```text
a = 0:  k = 1,
a = 1:  k = (1+eta^2)^(-2).
```

The `a=1` interval value uses exact rational reciprocal-square endpoint
monotonicity.

## Eta derivative recurrence

Differentiating

```text
(1+eta^2) k' + 4 a eta k = 0
```

`n` times gives, for full derivatives `d_n=k^(n)`,

```text
(1+eta^2) d_(n+1)
 + (2n+4a) eta d_n
 + (n(n-1)+4an) d_(n-1) = 0.
```

For normalized Taylor coefficients `c_n=d_n/n!` this becomes

```text
c_0 = k,
c_1 = -4 a eta c_0 / (1+eta^2),

c_(n+1) =
  -((2n+4a) eta c_n + (n-1+4a)c_(n-1))
   / ((n+1)(1+eta^2)),   n>=1.
```

The denominator is strictly positive on the real domain. Production uses this
identity with exact rational interval operations; it does not finite-difference
or fit sampled kernel values. `PressureKernelJetEnclosure` exposes the
normalized coefficients and the factorial-scaled full derivatives separately.

## Shared parameter boxes and subdivision

An exponent interval is not replaced by its midpoint. Every derivative row for
one sub-box reuses the same exponent and eta interval objects. Optional
rectangular subdivision evaluates rigorous sub-box enclosures and hulls them.
Thus refinement may reduce dependency overestimation without discarding the
shared parameter uncertainty. A finite `max_cells` cap prevents accidental
unbounded work.

## Independent checks recorded for this increment

The isolated focused command

```text
python -m pytest -q -W error tests/test_schedule_axis_pressure_kernel_enclosure.py
```

returned `15 passed in 0.21s`.

The source and focused test also passed

```text
python -m py_compile \
  openai_ns_reconstruction/schedule_axis_pressure_kernel_enclosure.py \
  tests/test_schedule_axis_pressure_kernel_enclosure.py
```

with exit code 0.

The regression checks exact `a=0` constant jets, exact `a=1` value/first/second
derivative identities at nonzero eta, factorial normalization, interval
parameter refinement/endpoint containment, an eta interval crossing zero,
entrypoint consistency, and fail-closed type/order/cell/range-reduction caps.

No full repository pytest, CLI demo, Lean build, or paper-exact audit result is
attributed to this increment in the local execution environment. Those checks
are recorded as `not run` unless repository CI later supplies exact-head
evidence.

## Real/complex boundary and remaining work

This constructor certifies only the real kernel and arbitrary requested finite
real eta derivative order. It does **not** certify a complex eta strip,
holomorphic extension, or a uniform all-order analytic norm. It also does not
enclose `clockWeight`, integrate the kernel over the full y-line, justify
differentiation under that integral, or certify the resulting axis pressure.
Those are downstream NS004/NS005/NS007/NS008 obligations.

`paper_exact_velocity_available=false` and `full_reconstruction=false` remain
mandatory.
