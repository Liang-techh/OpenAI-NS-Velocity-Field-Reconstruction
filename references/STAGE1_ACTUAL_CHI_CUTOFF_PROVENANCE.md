# Stage-1 actual-schedule chi-cutoff certificate provenance

## Scope

This increment is built directly from repository `main`
`3cc438f29524acb068e40bf3a99a8293c9222628` and stays strictly downstream of
the actual SchedulePressure low-|Z| margin. It does not implement
`AxisCoefficientSpace`, coefficient products, `naturalRemainder`, `x2`, a
fixed point, or materialized `phi/u/average/pressure`; Agent 6 retains that
backend lane.

The landed `schedule_axis_margin.certify_schedule_low_Z_margin(...)` already
constructs an analytic, no-sampling witness `m>0` such that on the actual
SchedulePressure datum,

`|Z(h,j,P,eta)| <= j/10  ->  m <= H(h,j,eta)^2`.

The landed `NaturalAxisRange` adapter then uses the pinned formal choice
`sigma=sqrt(m)/20`. The remaining downstream seam is exact arithmetic:
`sigma^2=m/400`, hence

`chi = H^2/(H^2+sigma^2) >= 400/401 > 99/100`.

`stage1_actual_chi_cutoff.py` binds those existing actual-schedule owners to
the exact rational constants `1/400`, `400/401`, `99/100`, and strict slack
`301/40100`. No eta samples, fitted margin, binary64 `0.99`, or caller-supplied
sigma/delta/chi threshold can substitute for the theorem identity.

## Pinned formal source

Official formal revision:
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant source:
`NavierStokes/NaturalAxisRange.lean`, especially
`low_Z_has_H_margin`, `exists_sigma`, and `exists_cutoff_parameters`.
The pinned theorem proves the strict implication

`|Z| <= j/10 -> 99/100 < chi`

after choosing `sigma = sqrt(m)/20`.

## Verification boundary

The regression uses the same theorem-admissible actual schedule fixture as the
landed margin and Picard paths (`P=2, m=1, lam=0.05, wait=30, h=0.01`,
`j=0.05`). It checks the actual analytic margin is consumed, the exact rational
chi algebra and strict slack are preserved, inexact downstream bounds are
rejected, and the constructor fails closed outside the actual margin
hypotheses.

The certificate does not turn the Python float realization of the analytic
margin into an interval or Lean proof object. It only removes the remaining
float-only chi-threshold metadata from the downstream handoff.

## Truth boundary

Status remains `formal-structure`. `full_reconstruction=false` and
`paper_exact_velocity_available=false` remain mandatory.

Still required: Agent 6's genuine complete `naturalRemainder(x1)`, `x2`,
global coefficient-space closure and fixed point; backend-owned norm-error
evidence; materialized `phi/u/average/pressure`; `NaturalProfileAssembly`;
regular inner-core / heat-exterior matching; and support/moment/matching/cone
closure.
