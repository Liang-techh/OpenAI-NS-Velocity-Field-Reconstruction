# Axis coefficient natural operator provenance

Status: **formal-structure**. `paper_exact_velocity_available=false` remains mandatory.

## Source pin

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Primary source: `NavierStokes/AxisResolvent.lean`
- Pinned definition: `AxisResolvent.naturalOperator`
- Supporting identity: `AxisReference.reference_coefficient`

The official operator is

`naturalOperator A = (1/2) * J_2 (chi * A)`,

with `J_2` the zero-datum regular radial inverse and `chi` the actual natural-axis multiplier.

## Executable construction

`src/openai_ns_reconstruction/axis_coefficient_natural_operator.py` builds the operator only from an `ActualScheduleReferenceAxisState`.

The actual `chi` parameter jets are not fitted or caller supplied. The pinned reference coefficient identity gives

`phi0[1] = -chi/4`,

so every actual eta derivative of the radially constant multiplier is recovered as `chi^(m) = -4 * phi0[1]^(m)` from the already-landed analytic actual-SchedulePressure angular reference family. Positive radial rows of `chi` are exactly zero.

The operator application then uses only the previously landed pinned coefficient primitives:

1. `coefficientOperators.product(chi, A)`;
2. `coefficientOperators.j2(...)`;
3. exact scalar multiplication by `1/2`.

No caller-supplied `sigma`, `epsilon`, radial divisor, coefficient table, derivative table, or replacement operator is accepted.

## Independent regression

`tests/test_axis_coefficient_natural_operator.py` checks:

- zeroth `chi` values against the actual schedule-derived `NaturalAxisData.chi` evaluator;
- the first eta derivative against an independent analytic derivative of `H^2/(H^2+sigma^2)`;
- exact radial constancy of `chi`;
- the coefficient identity `Q(A)[n] = (1/2) chi A[n-1]/(n(n+1))` and its first eta derivative on the actual angular reference state;
- the pinned resolvent equation on the already-landed closed-form reference coefficients, `phi0 + Q(phi0) = one` coefficient by coefficient;
- fail-closed rejection of a coefficient state on a different epsilon scale.

## Remaining boundary

This increment materializes only the pinned linear operator `Q`. It does **not** provide the complete weighted `AxisSpace` norm certificate, the infinite alternating series `naturalResolvent = sum_k (-Q)^k`, or a truncation/error certificate for evaluating that series on arbitrary coefficient states. Therefore the angular resolvent action required by `naturalRemainder`, `naturalRemainder(x0)`, the first genuine Picard iterate, fixed-point `phi/u`, derived average/pressure, `NaturalProfileAssembly`, and final support/moment/matching/cone checks remain unresolved.
