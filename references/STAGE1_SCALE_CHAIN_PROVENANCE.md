# Stage 1 actual-schedule scale-chain provenance

Status: **formal-structure / diagnostic certificate only**.  This file does not
promote the leading profile to paper-exact and does not claim that the true
natural resolvent norm is large.

## What is connected

`src/openai_ns_reconstruction/stage1_scale_chain.py` connects already-landed
actual SchedulePressure data to the existing coefficient-space bookkeeping:

1. `certify_actual_schedule_analytic_inputs(data, j)` supplies the theorem-side
   schedule `sigma`, an explicit common complex neighborhood `rho`, the common
   eleven-field value bound `B`, and a conservative compact-set
   `realPartSup(axisPhase)` bound.
2. The canonical `NaturalAxisCoefficients` choice gives `epsilon=rho/2`,
   `radiusLoss(1/2)=12`, hence the common coefficient-family norm upper bound
   `12 B`.
3. The pinned `AxisResolvent` filtration estimate uses
   `K = 2560 ||chi||`; with `||chi|| <= 12 B` the landed conservative ledger
   therefore uses an upward-rounded `K <= 2560*(12 B)`.
4. Before converting the complete factorial-series majorant to binary64, the
   new adapter checks positive terms
   `a_k = K^k/(k!(k+1)!)` using downward-rounded Decimal recurrence.

The official pinned Lean source is
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`,
`NavierStokes/AxisResolvent.lean`:

- `factorialMajorant (K) (k) = K^k/(k! (k+1)!)`;
- `factorialMajorant_nonneg` proves every term is nonnegative for `K>=0`;
- `factorialMajorant_succ` gives the recurrence
  `a_(k+1) = K/((k+1)(k+2)) * a_k`;
- `norm_alternatingResolvent_le` bounds the true alternating resolvent norm by
  the sum of these positive majorant terms.

## New verified boundary

For the existing concrete theorem-admissible regression schedule
`P=2, m=1, lambda=0.05, wait=30, h=0.01, j=0.05`, the actual-schedule analytic
certificate is now fed into this chain without caller-supplied `rho`, `B`,
resolvent norm, remainder constants, `Lambda`, or `C`.  A finite positive term
of the resulting **conservative majorant series** already exceeds
`sys.float_info.max`; the adapter records that term and stops before the
binary64 remainder/Lambda-C machinery.

This result is useful because it identifies the first concrete downstream
obstruction after PR #52: the present conservative schedule -> common-field ->
resolvent-majorant ledger is too large for the current numeric representation.
It is no longer merely unknown whether an overflow appears later in the chain.

## What this does *not* prove

The overflow witness is for the *upper majorant evaluated at the conservative
upper K*.  It is **not** a lower bound on the actual operator resolvent, is not a
failure of the Lean existence theorem, and is not evidence that the
coefficient-space fixed point does not exist.  In particular, replacing this
with a toy/sampled smaller K would be invalid.

The next theorem-faithful options are therefore to tighten one or more analytic
upper inequalities (especially the theorem-side sigma/common-field ledger), or
to carry the majorant/remainder selection in a representation that does not
require binary64, while preserving every contraction hypothesis.  Until that
is done, `phi/u/average/pressure`, `NaturalProfileAssembly`, and the final
support/moment/matching/cone checks remain unresolved and
`paper_exact_velocity_available` must remain false.
