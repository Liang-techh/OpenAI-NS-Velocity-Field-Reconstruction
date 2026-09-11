# Section 9 -> Section 10 endpoint-bridge provenance

Status: **formal-structure only**.

This increment does not construct the genuine Section 9 correction sequence,
does not prove an endpoint limit, and does not change
`paper_exact_velocity_available=false`.

## Sources and existing repository inputs

The Section 9 input is OpenAI, *Finite Time Blowup for Navier–Stokes*, Section
9.4, Lemma 9.8, especially Eqs. (9.17)--(9.19).  The repository's exact
rational exponent ledger is in `section9_residual_decay.py`, and PR #99 added
`section9_stage_certificate.py`, which checks separately certified norm data
against the paper envelopes at **one supplied q**.

The Section 10 endpoint input is the pinned formalization used by
`endpoint_limit_majorant.py` and `time_localization.py`.  A sufficient condition
for a locally uniform endpoint limit of `D^n R` on one compact spatial window is
an integrable bound

`sup_x ||partial_t D^n R(t,x)|| <= C (1-t)^(-alpha)`,

with `0 <= alpha < 1`, on a late interval ending at `t=1`.  The official time
switch is exactly one on the executable representative for `t >= 3/4`, so a
majorant whose full validity interval starts there transfers through time
localization without cutoff-correction terms.

Canonical source URLs and the pinned OpenAI snapshot are recorded in
`references/SOURCES.md` and the existing provenance files.

## What this increment adds

`section9_endpoint_bridge.py` records the first explicit fail-closed interface
between those two existing layers.

`Section9UniformEndpointDerivativeWitness` requires a separately justified
**uniform** estimate on the whole stated spatial compact and whole late-time
interval.  For an endpoint jet of total spacetime degree `n`, it requires the
Section 9 residual derivative order to be exactly `n+1`, because
`partial_t D^n R` is one order higher.  This derivative-order check prevents an
off-by-one use of Lemma 9.8.

The witness also requires:

- a nonnegative coefficient with an evidence classification and nonempty
  provenance via `CertifiedBoundDatum`;
- an integrable exponent `0 <= alpha < 1`;
- the official Section 10 endpoint `t=1`;
- a validity interval beginning no earlier than the exact `t>=3/4` unit
  plateau; and
- an explicit Section 9 stage and compact spatial-window identifier for audit.

`admit_section9_uniform_endpoint_witness` turns only such a witness into the
existing `EndpointPowerLawMajorant` and immediately passes it through
`section10_endpoint_localization_transfer`.

A `Section9PointwiseBoundCertificate` is rejected explicitly.  Passing Eq.
(9.18) at one `q` is not a proof of a bound uniform in `t` and `x`.

## Why no automatic q -> (1-t) conversion is made

The similarity chart satisfies

`1-t = q (1-eta^2)`

and, equivalently in physical variables,

`q - z^2 q^(2h) = 1-t`.

On a fixed physical compact this does **not** give a global upper comparison
`q <= C (1-t)^gamma`: away from the axis, `q` need not approach zero with
`1-t`.  Therefore a pointwise `q`-power from Eq. (9.18) cannot safely be
rewritten as an endpoint time majorant without the genuine Section 9 support /
common-domain theorem.  This module intentionally refuses to invent that
missing comparison.

## Independent regression checks

`tests/test_section9_endpoint_bridge.py` verifies four distinct boundaries:

1. an order-`n+1` uniform witness enters the official Section 10 late plateau,
   and its endpoint-tail budget is checked against the closed antiderivative
   `C/(1-alpha) (1-t)^(1-alpha)`;
2. an off-by-one Section 9 derivative order and a transition-window start fail
   closed;
3. nonintegrable exponents and a non-`t=1` endpoint fail closed; and
4. a genuinely passing PR-#99 Eq. (9.18) pointwise certificate is still rejected
   when presented as an endpoint-majorant witness.

The fourth regression is deliberately negative: it protects the repository
from turning a valid local arithmetic check into a false uniform theorem.

## Remaining blocker

This bridge becomes substantive only when the actual Section 7--9 construction
provides theorem-derived or rigorous interval witnesses for
`sup_x ||partial_t D^n R||` on the common late-time domain for every `n`.
Equivalently, an upstream proof may provide the missing support/common-domain
comparison that converts Proposition 9.9 flatness into these integrable
endpoint majorants.

Until then the record keeps
`source_theorem_machine_verified=false`,
`actual_section9_sequence_verified=false`,
`endpoint_limit_constructed=false`, and
`paper_exact_velocity_available=false`.
