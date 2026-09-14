# Finite physical-power ladder provenance

Scope: Agent 7 / Issue #2, downstream of the hierarchy-bound recursive SlowBorel cutoff, exact retained-recurrence cancellation, and coherent finite residual-extension gates in PR #290.

## Source chain

The certificate consumes only already-certified `FiniteCoherentSlowBorelResidualChain` endpoints. Each endpoint carries a `CertifiedFullResidualTailMajorant` whose first-omitted exponent is the physical slow-weight exponent implemented from the pinned Section 5 residual algebra:

`base_power + 2 * (N + 1) * h`.

The implementation rechecks that identity using exact `Fraction` arithmetic (`slow_order_exact`) at every endpoint and rechecks the exact `2*h` exponent gain between adjacent truncation orders. It does not introduce a generic cutoff, new coefficient recurrence, sampled cancellation, or caller-supplied `C[j,m]` row.

## New admission rule

A finite tuple of requested q-powers is admitted only when:

- it contains exactly one target per coherent truncation endpoint;
- every target is an exact positive integer/Fraction (floats are rejected);
- targets are strictly increasing;
- the regime has `0 < q < 1`;
- every stored first-omitted exponent exactly matches the physical slow-weight formula;
- each requested target is no larger than the certified first-omitted exponent; and
- the omitted-tail coefficient prefactor and residual log-majorant remain finite/non-NaN.

For each endpoint the certificate also exposes a normalized log-majorant for `residual / q**target`, so coefficient-prefactor size is not discarded when a target exponent is reported.

## Truth boundary

This is a finite-prefix target-power certificate only. It explicitly reports:

- `finite_prefix_only = true`
- `infinite_coherent_family = false`
- `uniform_in_q = false`
- `all_jets_flat = false`
- `super_algebraic = false`
- `paper_exact = false`

It therefore does **not** prove Proposition 5.3 or arbitrary-order/super-algebraic convergence. The remaining upstream obstruction is the same one exposed by PR #290: Agent 2 must provide, for every requested hierarchy order, hierarchy-owned analytic derivative/support bounds and exact retained recurrence row decompositions under one uniform all-order envelope. Agent 2 PR #301 currently advances the real mixed-jet hierarchy but still fail-closes beyond the finite implemented derivative frontier, so this downstream certificate does not infer missing higher-order rows.

## Regression coverage

`tests/test_background_target_decay.py` checks:

1. exact reconstruction of each physical first-omitted exponent and retained coefficient-prefactor normalization;
2. rejection of a requested power larger than the certified first-omitted exponent;
3. rejection of floating or non-increasing target powers; and
4. rejection of `q = 1` as a non-decaying regime.
