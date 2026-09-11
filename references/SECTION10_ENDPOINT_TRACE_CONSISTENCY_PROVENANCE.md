# Section 10 endpoint-trace consistency provenance

Status: **formal-structure only**.

This increment adds an independent finite-order falsification gate for the
conditional endpoint-majorant path. It does **not** derive endpoint estimates
from Section 9, does not construct endpoint limits, and does not change
`paper_exact_velocity_available=false`.

## Source and landed interfaces

The pinned OpenAI formalization's `CandidateFromLimits` requires locally uniform
limits of every full spacetime derivative of the actual closed-past
Navier--Stokes residual as `t -> 1-`. The landed
`endpoint_limit_majorant.py` records one sufficient, caller-supplied hypothesis,

`sup_x ||partial_t D^n R(t,x)|| <= C (1-t)^(-alpha)`, with `0 <= alpha < 1`,

and integrates it to a uniform Cauchy budget. The landed
`section10_endpoint_ladder.py` packages degrees `0..N` on the official
`t >= 3/4` time-switch plateau, while `section10_endpoint_trace.py` uses the
same hypotheses to enclose hypothetical endpoint jets.

## What this increment adds

`section10_endpoint_trace_consistency.py` checks a logically different necessary
consequence. At one spatial point and two certified pre-endpoint times it
evaluates the supplied dense full jets, with shape `(3,) + (4,)*n`, and requires

`||D^n R(t1,x) - D^n R(t0,x)|| <= integral_{t0}^{t1} C(1-s)^(-alpha) ds`

for every degree in the finite ladder.

The right-hand side comes only from the independently supplied analytic
majorant. The implementation never estimates, fits, or shrinks `C` or `alpha`
from the sampled jets. A violation raises and therefore cannot be relabeled as a
successful endpoint certificate.

The check is pointwise and finite-order. Passing it is **not** evidence that the
majorant holds uniformly on the spatial compact, nor that the queried jets are
the actual Section 9 residual jets.

## Independent regression

`tests/test_section10_endpoint_trace_consistency.py` uses the analytic family

`J_n(t) = L_n + c_n sqrt(1-t) e_n`.

Its derivative norm is exactly `(c_n/2)(1-t)^(-1/2)`, so the production
integrated budget can be compared against the independently evaluated closed
form

`c_n (sqrt(1-t0) - sqrt(1-t1))`.

A second test injects an extra component at the later sample. The verifier must
reject the data rather than refit the majorant. Shape mismatches and times
outside the common late-time window also fail closed.

## Boundary retained

The following remain unverified and are deliberately represented as false in
the returned aggregate certificate:

- that the supplied jet family is the actual localized Section 9 residual;
- that the power-law majorants hold uniformly on the indexed spatial compact;
- locally uniform endpoint limits;
- arbitrary-order endpoint control and smooth Borel gluing;
- smooth compact forcing, bounded kinetic energy, and blow-up closure;
- `paper_exact_velocity_available`.

The next substantive upgrade still requires machine-connected Section 9
correction/residual values together with a proof of uniform physical-time
derivative majorants up to `t=1`.
