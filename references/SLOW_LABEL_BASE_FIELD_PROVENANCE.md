# Section 6 slow-label/base-field bridge provenance

Status: **formal-structure**. This file does not upgrade Stage 3 to paper-exact.

## Paper mapping

The implementation in `src/openai_ns_reconstruction/slow_labels.py` is mapped to the published construction as follows.

- Section 6.2, Eq. (6.8): labels `gamma=(ell,a,sigma)`, the active dyadic shell `1/2 <= q/Q <= 2`, and `X_a <= r^2/(2q) <= X_b`.
- Section 6.2, Eq. (6.9): the two signs share the same slow box/cutoff. The runtime object only certifies a mesh-scale support enclosure; it does **not** manufacture the missing squared partition `chi_ell chi_{ell,a}`.
- Section 6.2, surrounding Eq. (6.8): the slow product mesh is `S_*^-3` in `(R,Z,T)` and representatives are fixed before derivatives.
- Eq. (6.11): once the Lemma 6.1 rectangle radius `r0` is supplied, `L_s=2 r0/c_i`.
- Section 7.1 and Eq. (7.2): the provider-supplied frozen base jet is converted to `F0`, `g0=(R0 F_R0,G_R0)`, the positive-growth frame, `k=ceil(epsilon^-1/2)`, and `B_s^2=lambda0/[epsilon k^2(1+u_*^2)^(3/2)]`.

The normalized active-shell calculation is not a fitted surrogate. Writing `q=Q s`, `z=Q^(1/2-h) Z`, and `tau=Q T` in the similarity identity gives `s-Z^2 s^(2h)=T`, so `q/Q=s` and `X=R^2/(2s)`. The code reuses the repository's robust scalar similarity solve for this dimensionless equation.

## Interface boundary

`TangentialBaseJetProvider` is intentionally only a typed hook. A future paper-exact Proposition 5.5 background can implement it and supply the normalized `F,G` jet at each chosen representative. No `paper_exact=True` flag exists and no caller-provided jet is promoted by this module.

`SlowBoxEnclosure` checks only the paper's `S_*^-3` support-size envelope and that the representative lies inside it. It is not evidence that the actual normalized squared partition of unity has been constructed. Likewise `rectangle_radius` and `u_star` have no arbitrary defaults: they must come from the upstream Lemma 6.1 / Section 7 choices.

## Independent tests

`tests/test_slow_labels.py` reconstructs a manufactured dimensionless active-shell point from a prescribed `q/Q`, checks the derived `X`, verifies failure outside `[1/2,2]`, checks sign duplication and the exact mesh exponent, and verifies that label phase data are frozen from exactly one provider call at the representative. The expected growth eigenvalue, pulse length, carrier integer, and `B_s` are recomputed independently in the test.

## Still missing for paper-exact status

- the actual smooth squared partitions `chi_ell` and `chi_{ell,a}` and their support certificates;
- the Lemma 6.1 common rectangle radius/centers and cross-label support separation;
- a paper-exact Proposition 5.5 background implementation of `TangentialBaseJetProvider`;
- analytic `LocalBaseBounds` / C1-C2 certificates on every enlarged active slow box, rather than point samples;
- the resulting uniform Eqs. (7.9)-(7.11) estimates, stress-cone realization, amplitude/curl waves, mean corrections, and Section 9 iteration.
