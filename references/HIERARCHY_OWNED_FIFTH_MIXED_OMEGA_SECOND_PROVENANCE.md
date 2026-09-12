# Section 5 hierarchy-owned fifth-mixed U / Omega second-eta provenance

## Scope

This increment closes the ownership gap left after the compact Lemma 5.2
fifth-mixed repair landed. It does **not** introduce a new Section 5 formula.

`Section5FifthMixedCoefficientJetSource` owns one authoritative repaired
`AxialFifthMixedJet`. Its fourth- and third-mixed U providers are exact
projections of that same object. `Section5LowerHistoryFifthMixedHierarchy`
checks the projection before exposing the lower layer and derives
`RegularFluxFourthMixedJet` only through the already-landed analytic Eq. (5.2)
mapper.

`hierarchy_owned_omega_second_parameter_jet` then feeds hierarchy-derived
fourth-mixed beta jets and hierarchy-projected axial second jets into the
already-landed analytic second-eta Eq. (5.6) implementation. It accepts no
caller beta/Omega derivative table.

## Paper / implementation map

- Lemma 5.2 and Appendix A: compact five-bump moment repair.
- Eq. (5.2): regular radial flux `beta_n = V_n / X`, including the landed
  positive-order `D + lambda_n` correction.
- Eqs. (5.3)-(5.6): strict-lower recursive source structure and the regular
  `Omega_k / X` row.
- Eq. (5.7): downstream PositiveAxis solve that will consume the resulting
  stronger forcing jet.

The new ownership layer reuses:
- `background_moment_repair_fifth_mixed_jets.py`;
- `background_regular_flux_fourth_mixed_jets.py`; and
- `background_omega_second_parameter_jet.py`.

It does not duplicate their formulas.

## Verification

Regression uses the real `Lemma52MomentRepair.from_intervals(...)` five-bump
repair with analytic degree-five moment jets.

1. Repaired order-1 hierarchy-owned fifth-mixed U is compared exactly with an
   independently instantiated `Lemma52RepairedFifthMixedJetAdapter`.
2. Its fourth/third projections are required to equal the inherited hierarchy
   paths exactly.
3. Hierarchy-owned beta fourth-mixed data are compared exactly with the direct
   landed Eq. (5.2) adapter.
4. The hierarchy-owned second eta derivative of `Omega_1 / X` is compared
   against a centered eta difference of the previously landed
   hierarchy-owned analytic first-eta path, both at `X=0` and inside an active
   compact repair support. Finite difference exists only as a test oracle.
5. A missing leading fifth-mixed U provider fails closed, and an inconsistent
   fifth-to-fourth projection is rejected.

## Truth boundary

Status remains **Stage 2 / `formal-structure`**.

The unrepaired fifth-mixed leading profile is still an upstream contract.
Issue #1 must supply the genuine order-zero strong profile before this bridge
can be used as a paper-backed recursive solve. Test fixtures, caller-supplied
base jets, generic cutoffs, and finite-order numerical fits are not
paper-exact.

This increment does not claim Picard convergence, final positive-order
coefficient materialization, a recursive cutoff-scale schedule, Proposition
5.3 all-jets residual decay, or a paper-exact reconstructed velocity.
