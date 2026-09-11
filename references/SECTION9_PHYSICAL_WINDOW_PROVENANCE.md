# Section 9 physical-slab local-sum provenance

Status: **formal-structure**. This note records one exact bridge between the landed Section 9 Eq. (9.21) local-finiteness admission and the fixed Section 10 spatial support. It does not construct the correction sequence, prove Proposition 9.9, or provide an endpoint residual limit.

## Paper / formalization map

- Similarity geometry: Eq. (4.1), equivalently `tau = q - z^2 q^(2h)` with `tau=1-t` and `0<h<1/2`.
- Shrinking-cutoff/local-finiteness mechanism: Lemma 5.4.
- Common small-q domain and local sum: Lemma 9.7, Proposition 9.9 Step 2, Eq. (9.21).
- Final fixed spatial localization: Section 10 / the pinned `SpatialLocalization.supportCylinder`, with `|z| <= 1/4` and hence `z^2 <= 1/16`.
- Final time localization late plateau: the pinned Section 10 switch is exactly one for `t >= 3/4` in the executable representative already on main.

## Exact physical-window bound

Let `t0 <= t <= t1 < 1` with `t0 >= 3/4`, and let the point lie in the fixed Section 10 support cylinder. Then `tau=1-t` satisfies `0 < 1-t1 <= tau <= 1-t0 <= 1/4`, while `z^2 <= 1/16`.

The physical similarity root obeys

`tau = q - z^2 q^(2h)`.

If `q >= 1`, then because `0<2h<1`, `q^(2h) <= q`; therefore `tau >= q(1-z^2) >= 15/16`, contradicting `tau <= 1/4`. Hence `q<1`, so `q^(2h) <= 1`, and consequently

`1-t1 <= q <= 1-t0 + 1/16`.

The implementation carries the time bounds and q-strip endpoints as exact `Fraction` values. It invokes `certify_eq_9_21_local_finiteness` only when the analytic upper endpoint is **strictly** below the independently supplied theorem-side `q_big`; no sampled q maximum is allowed to enlarge the domain.

## Independent regression

`tests/test_section9_physical_window.py` separately evaluates the existing Eq. (4.1) physical-root solver over several `h`, `t`, and `z` values spanning the official support/late slab and checks that every computed q lies inside the analytic exact-rational enclosure. Additional tests verify that:

- moving the physical slab later tightens the analytic q upper bound;
- a theorem-side `q_big` that is too small is rejected rather than inferred from samples;
- `t1=1` is rejected because the uniform positive q lower bound collapses;
- slabs beginning before the official `t>=3/4` plateau are rejected.

The numerical grid is only an independent regression of the algebraic bound; it is not the source of the certificate.

## Boundary deliberately preserved

This bridge is finite-slab infrastructure only. It does **not** show that one fixed positive q-strip covers the singular endpoint: at the origin `z=0`, Eq. (4.1) gives `q=1-t`, so `q -> 0` as `t -> 1-`. Therefore the certificate cannot be converted into the uniform physical-time derivative majorants required by the Section 10 endpoint-limit/Borel chain.

Still missing are the genuine Section 7/8 correction fields, the paper-derived infinite scale schedule tied to those fields, the constructed Eq. (9.21) sum, arbitrary-order Eq. (9.20) residual flatness, the locally uniform all-order endpoint limits, smooth force extension through `t=1`, and the final bounded-energy/blow-up closure. `paper_exact_velocity_available` remains `false`.
