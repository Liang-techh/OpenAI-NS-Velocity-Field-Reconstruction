# Field definition and fixed contract

The executable definitions are in `src/ns_reconstruction/_runtime/`. Their arithmetic is imported from the pinned study implementation; only sibling imports are namespaced. MATLAB evaluates the same basis after baking in the QR transforms and initial-energy normalization.

## Axisymmetric representation

Set `s = x^2 + y^2 = r^2`. A poloidal potential `F(s,z,t)`, swirl factor `B(s,z,t)` and pressure `p(s,z,t)` define

```text
A = -F_z
C = 2 F + 2 s F_s
u = (x A - y B, y A + x B, C)
```

Consequently `2 A + 2 s A_s + C_z = 0` algebraically. This divergence identity is not the momentum equation. The streamfunction is `psi = r^2 F`; `u_theta = r B`. Opposite axial-parity basis columns allow nonzero midplane axial bias/shear without breaking axisymmetry.

Each block has `9 * 12 * 8 = 864` coefficients. There are 1,728 velocity coefficients, 864 pressure coefficients and two force parameters, for 2,594 stored parameters. The spatial basis uses smooth bump-weighted Legendre columns, derivative-enriched columns and opposite-parity axial columns. Time is represented by eight Chebyshev terms `T_k(4t - 2)`. No time interpolation is used by the spectral evaluators.

The bump is `beta(q) = exp(1 - 1/(1-q))` for `q < 1` and zero otherwise. The spatial bump factors use `q=s/4` and `q=z^2/4`. Their mathematical extension is smooth at the support boundary; very small floating-point tails are handled by the pinned implementation. This does not establish every derivative by interval arithmetic.

## Prescribed force

Let `b(s,z) = beta(s/4) beta(z^2/4)` and `g(t) = beta((2t-1)^2)`. With constants `a,c` in `[0,10]`, the independently specified force is

```text
A_f = -a (b + z b_z)
B_f = c (b + s b_s)
C_f = 2 a z (b + s b_s)
f = g(t) (x A_f - y B_f, y A_f + x B_f, C_f)
```

This family is divergence-free. Its parameters are frozen with the candidate; the force is not created after the fact from its residual.

## Complete momentum residual

```text
R = u_t + (u dot grad) u + grad p - 0.01 Delta u - f
```

The analytic evaluator includes all three Cartesian components, convective terms, pressure gradient, viscosity and the prescribed force. Independent acceptance recomputes derivatives from field values using Cartesian finite differences.

## What these fields do not establish

The publication window ends at `t=0.75`. The package does not extrapolate to a singular endpoint or prove unbounded vorticity, global regularity, existence of an exact nearby solution, or correspondence with an unavailable original velocity field. A finite polynomial time representation and favorable plots are not a blow-up construction. ST054 changed pressure relative to its parents, not their velocity geometry. Both complete momentum gates still fail.
