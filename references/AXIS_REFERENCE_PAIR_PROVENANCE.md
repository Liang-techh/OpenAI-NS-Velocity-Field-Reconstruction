# Stage-1 natural-axis reference-pair coefficient provenance

Status: **formal-structure**. This increment materializes genuine coefficient functions for the pinned reference pair, but it does **not** materialize the nonlinear Theorem 4.6 fixed point and does not change `paper_exact_velocity_available=false`.

## Pinned official source

This increment is mapped to `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`:

- `NavierStokes/AxisContraction.lean`
  - `referencePair`
  - `coefficientOperators`
- `NavierStokes/AxisReference.lean`
  - `reference_coefficient_zero`
  - `reference_coefficient_succ`
  - `reference_coefficient`
- `NavierStokes/AxisWeightEstimates.lean`
  - `radialDivisor`
  - `regularInverseJet`
- `NavierStokes/NaturalAxisCoefficients.lean`
  - `window`
  - fixed fields `inverseL`, `zStar`, and `chi`
- `NavierStokes/NaturalAxisData.lean`
  - `Z`
- `NavierStokes/SchedulePressure.lean`
  - actual axis pressure datum used inside `zStar`

The pinned contraction uses

`x0 = (S one, -(1/2) J_1(inverseL * zStar))`,

where `S` is the natural angular resolvent. `AxisReference.reference_coefficient`
proves the angular radial coefficients explicitly:

`phi0[n](eta) = (-chi(eta)/2)^n / (n! (n+1)!)`.

For the axial component, `inverseL*zStar` is radially constant and
`regularInverseJet 1` has coefficient zero at radial degree zero and
`source[n]/radialDivisor(1,n)` at degree `n+1`. Hence the reference axial
component has only

`u0[1](eta) = -(1/2) inverseL(eta) zStar(eta)`

nonzero, because `radialDivisor(1,0)=1`.

## Landed executable capability

`src/openai_ns_reconstruction/axis_reference_pair.py` now constructs these
reference coefficients from the **actual landed SchedulePressure chain**:

1. `certify_schedule_low_Z_margin(data,j)` supplies the theorem-side
   `sigma=sqrt(margin)/20`; there is no independent caller `sigma` on the
   factory path.
2. `chi(h,j,sigma,eta)` supplies the radially constant coefficient used by the
   actual natural resolvent.
3. The actual `SchedulePressure` pressure and its proved derivative formula are
   inserted into `NaturalAxisData.Z` on the full coefficient window
   `[-11/10,11/10]`.
4. The angular and axial reference coefficients are then evaluated by the
   pinned resolvent/regular-inverse recurrences above.

This closes a concrete part of the previous “referencePair unavailable”
blocker: the zeroth parameter-jet radial coefficient functions of both members
of `x0` are now executable and tied to the actual schedule rather than an
arbitrary coefficient array.

## Independent checks

`tests/test_axis_reference_pair.py` checks the angular recurrence implementation
against the separate factorial closed form proved in `AxisReference`; verifies
the resolvent recurrence coefficient-by-coefficient; rebuilds
`NaturalAxisData.Z` independently from actual schedule pressure and checks the
axial reference against the regular integrated equation
`2*(Y u'' + u') = -inverseL*zStar`; and confirms that the evaluator covers the
full pinned `NaturalAxisCoefficients.window`, not just `[-1,1]`.

## Explicit non-claims / remaining blocker

This increment is **not** a Python implementation of the complete Banach space
`AxisCoefficientSpace.AxisSpace`. It materializes coefficient functions
(`m=0`) only; it does not yet store/certify every compatible parameter jet, and
therefore cannot yet instantiate the generic bounded operators as a complete
normed-space backend.

It also does not evaluate `naturalRemainder` or claim that one Picard step has
been executed. The next theorem-faithful step is to lift these actual reference
coefficient functions into a compatible coefficient-state representation with
analytic parameter jets, implement the finite coefficient actions needed by
`coefficientOperators`, and evaluate the first real `naturalRemainder(x0)`.
The already-landed `axis_fixed_point_picard.py` can then apply the certified
`1/(2 Lambda)` update and geometric tail budget.

The actual schedule pressure evaluator uses the repository's documented
analytic reductions plus finite Gauss-Legendre quadrature on smooth transition
intervals. Consequently the numerical `zStar` values here inherit that caveat;
this increment is not interval arithmetic and is not a Lean proof object.
