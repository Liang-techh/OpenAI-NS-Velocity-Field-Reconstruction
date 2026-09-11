# SchedulePressure low-|Z| margin provenance

Truth status: **formal-structure subconstruction**.  This file records an explicit
quantitative realization of the compactness step used by the pinned Lean theorem
`NaturalAxisRange.low_Z_has_H_margin`.  It is tied to the actual outgoing
`SchedulePressure.axisPressure`; it is not the earlier zero-positive-time
ideal-prefix witness and it does not complete Theorem 4.6.

## Upstream pin

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Relevant files:
  - `NavierStokes/NaturalAxisData.lean`
  - `NavierStokes/NaturalAxisRange.lean`
  - `NavierStokes/OutgoingSchedule.lean`
  - `NavierStokes/OutgoingTail.lean`
  - `NavierStokes/SchedulePressure.lean`

The formal chain used here is:

1. `SchedulePressure.pressureData` instantiates `NaturalAxisData.PressureData`
   for the actual schedule when the ideal-prefix amplitude satisfies `P >= 2`.
2. `NaturalAxisRange.exists_unique_root` gives the unique root
   `r in (-j/4,-j/5)` of `H`.
3. `NaturalAxisRange.Z_at_root_lower` gives `Z(r) > j/5`.
4. `NaturalAxisData.cross_identity` relates `H(eta)` to distance from `r`.
5. `NaturalAxisRange.exists_sigma` chooses `sigma=sqrt(m)/20` once a positive
   uniform low-|Z| margin `m` is supplied.

The Lean proof obtains `m` abstractly by compactness.  The runtime module
`src/openai_ns_reconstruction/schedule_axis_margin.py` replaces only that
nonconstructive minimum-selection step with the following conservative explicit
inequality chain.

## 1. Clock-mass envelope

Write `M0 = integral_R clockWeight(y) dy`.  The schedule has exact ideal prefix
mass `5 P^2`.  For `0 <= y <= releaseStart`, the release and terminal-tail
factors are inactive and `logAmplitude <= 1/10`; therefore

`clockWeight(y) <= exp(1/5) P^2 < (5/4) P^2`.

For `y >= releaseStart`, the combined `logAmplitude + releaseAdjustment`
derivative is at most `-(1/2+h)`.  This uses monotonicity of the pinned Gaussian
step and `2h < lam`.  The terminal taper contributes at most `1/(1-rho)`, with

`rho = exp(-5) h / (16 * 33) < h / 16896`,

because the repository has already fixed the theorem-admissible derivative-bound
witness `S=32` and `exp(-5)<1/32`.

The code therefore uses

`M = 5P^2 + (5/4)P^2 R + (5/4)P^2 / ((1-rho_bar)^2 (1+2h))`,

where `R` is an explicit upper bound for `releaseStart` and
`rho_bar=h/16896`.  The release bound uses only printed schedule parameters:

`releaseStart = exp(m)+13+wait+13/lam+330 log 2+30 log(1/lam)`.

It substitutes `log 2 < 1`, `log(1/lam)<1/lam`, and the elementary inequality
`exp(m) < (1-m/n)^(-n)` for integer `n>m`, with
`n=2 ceil(m)+2`.  No tailDebt, decayHold, numerical quadrature, or sampled
clock-weight maximum enters this bound.

## 2. Pressure derivative envelope

For the actual pressure kernel

`K_a(eta)=(1+eta^2)^(-2a)`, `0<=a<=1`,

and `|eta|<=1`, positivity gives

- `|P| <= M/2`,
- `|P'| <= M`,
- `|P''| <= 4M`.

The first derivative is the pinned `SchedulePressure.axisPressure_hasDerivAt`
formula.  The second inequality follows by differentiating its positive
exponent-weighted integral once more and bounding both exponent moments by the
total clock mass.

## 3. Global Z Lipschitz bound

Using

`Z = -A(1-2 eta U)U - 4H - d P' + 4A eta P`,

`U=4eta+j`, `A=1/2+h`, `D=1/2-h`, the implementation bounds the polynomial
part of `|Z'|` on `[-1,1]` by

`A [2(8+j)(4+j) + 4(9+2j)] + 4[D + 2(4+j) + 4]`.

The pressure-dependent part is at most `(6+6A)M`.  Their sum is the recorded
`z_lipschitz_upper = L_Z`.

## 4. Explicit H margin

At the unique H-root `r`, the formal theorem gives `Z(r)>j/5`.  Thus any
`eta` with `|Z(eta)|<=j/10` satisfies

`|eta-r| > j/(10 L_Z)`.

At the same root, the exact cross identity reduces to

`H(eta)d(r) = (eta-r)[D(1+eta r)+4d(eta)d(r)]`.

Because `r in (-j/4,-j/5)`, `|eta|<=1`, and `d>=0`, one obtains

`|H(eta)| > D(1-j/4)|eta-r|`.

The runtime chooses the deliberately smaller non-strict witness

`m_H = 1/2 * [D(1-j/4) * j/(10 L_Z)]^2 > 0`.

That value is then passed unchanged to the already-landed
`cutoff_parameters_from_margin`, yielding the pinned theorem choice
`delta=j/10`, `sigma=sqrt(m_H)/20` and `chi>99/100` on the low-|Z| set.

## Verification and boundary

`tests/test_schedule_axis_margin.py` independently checks that the analytic
clock-mass envelope dominates the executable all-real-line schedule mass, that
the actual schedule `Z` has the expected `j/5` gap at the numerically located
formal H-root, and that sampled low-|Z| points obey the constructed margin.
Those samples are regression checks only; they do not construct `m_H`.

The Python runtime evaluates the explicit real-arithmetic bound in binary64 and
is not an interval-arithmetic or Lean certificate.  Therefore this increment
closes the **missing explicit margin-selection formula** for a concrete
`TailData`, but does not promote Stage 1 to `paper-exact`.  Still missing are the
paper's constructive `Lambda/C` choices, coefficient-space fixed-point fields
`phi/u/average/pressure`, `NaturalProfileAssembly` instantiated from those real
fields, and the final support/moment/matching/cone certificates.
