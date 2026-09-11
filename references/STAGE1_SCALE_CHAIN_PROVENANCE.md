# Stage 1 actual-schedule scale-chain provenance

Status: **formal-structure / diagnostic certificate only**.  This file does not
promote the leading profile to paper-exact and does not claim that the true
natural resolvent norm is large.

## Pinned source

Official source commit:
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant pinned modules are
`NavierStokes/NaturalAxisCoefficients.lean` and
`NavierStokes/AxisResolvent.lean`.

`NaturalAxisCoefficients.boundedAxisElement_norm` is a one-field statement:
if a complex field is analytic on the common tube and has value bound `B_k`,
then at coefficient radius `epsilon<rho` its AxisSpace norm is bounded by

`B_k * radiusLoss(epsilon/rho)`.

For the canonical `epsilon=rho/2`, the already-landed exact identity is
`radiusLoss(1/2)=12`.  The pinned existence proof chooses one common `B` first
and therefore obtains the convenient common bound `12 B`, but the theorem does
not require the same value bound to be used for every call to
`boundedAxisElement`.  The resulting eleven elements still fit the same
`CoefficientFamily` structure by taking the maximum of their individual norm
bounds as the structure's common `bound` field.

`AxisResolvent.lean` then uses the norm of the **chi element specifically**:

- `factorialMajorant K k = K^k/(k! (k+1)!)`;
- `naturalOperator_pow_bound` uses `K = 2560 * ||chi||`;
- `norm_alternatingResolvent_le` bounds the resolvent by the sum of those
  positive factorial-majorant terms.

## Componentwise refinement

`src/openai_ns_reconstruction/axis_componentwise_input_bounds.py` applies the
same pinned one-field Cauchy estimate separately to the eleven value bounds
already produced by `schedule_analytic_neighborhood.py`:

`||element_k|| <= 12 B_k`.

It records

- the common coefficient radius `epsilon=rho/2`;
- every individual fixed-field norm upper bound `12 B_k`;
- the common `CoefficientFamily.bound` upper as `max_k 12 B_k`;
- the normalized amplitude bound `M=12`;
- the chi-specific norm upper `12 B_chi`;
- `K <= 2560 * (12 B_chi)` and the complete positive factorial-series
  enclosure.

This is not a sampled or hand-tuned reduction.  It removes only an avoidable
cross-field overestimate: a large bound for `complexGradient` no longer enters
`AxisResolvent`, because the pinned resolvent theorem depends on chi rather than
on the maximum norm of all fixed fields.

## Actual-schedule consequence

`src/openai_ns_reconstruction/stage1_scale_chain.py` now uses that componentwise
ledger for the actual landed SchedulePressure certificate.  On the existing
regression schedule

`P=2, m=1, lambda=0.05, wait=30, h=0.01, j=0.05`,

the chi-specific majorant parameter is below `100000`, and the complete
factorial resolvent majorant is representable in binary64.  Thus the previous
"one positive factorial-majorant term already exceeds binary64" obstruction was
an artifact of feeding the unrelated common eleven-field bound into chi; it is
no longer the first blocker on this actual schedule.

The chain now proceeds through the factorial-series enclosure and fails closed
at the next representability boundary during conservative
`remainderBound/remainderLip` propagation.  The diagnostic records that later
stage explicitly instead of manufacturing finite values.  This new overflow is
again only about the current conservative upper-bound arithmetic; it is not a
lower bound on the true remainder, not a failure of the Lean existence theorem,
and not evidence against existence of the fixed point.

The standalone downward-rounded positive-term witness remains in the module for
large-K cases: if such a term exceeds `sys.float_info.max`, the code still stops
before converting the series.

## Truth boundary

This increment does **not** materialize the coefficient-space Banach elements,
`phi`, `u`, `average`, or `pressure`; it does not construct
`NaturalProfileAssembly`; and it does not close support, moments, matching, or
cone conditions.  It only tightens a theorem-faithful analytic norm ledger and
moves the actual-schedule diagnostic to the next concrete blocker.

`paper_exact_velocity_available` therefore remains false and Stage 1 remains
`formal-structure`.
