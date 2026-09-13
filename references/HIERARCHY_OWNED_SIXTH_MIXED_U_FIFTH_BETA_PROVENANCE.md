# Section 5 hierarchy-owned sixth-mixed U / fifth-mixed beta provenance

## Scope

This increment closes one derivative-ownership seam on the route to the
hierarchy-owned third eta derivative of the strict-lower source. The landed
second-eta Eq. (5.6) path owns a fourth-mixed regular flux `beta_n = V_n / X`.
Differentiating the Eq. (5.6) row once more requires the three total-order-five
regular-flux entries `beta_XXetaetaeta`, `beta_Xetaetaetaeta`, and
`beta_etaetaetaetaeta`.

Eq. (5.2) derives those rows from one stronger axial layer:
`U_XXetaetaetaeta`, `U_Xetaetaetaetaeta`, and `U_etaetaetaetaetaeta`.
The implementation baseline is repository `main`
`1d812ff01d4575551f804c3599d6efa8985c7d2f`, which already includes PR #283.

## Analytic Eq. (5.2) lift

`background_regular_flux_fifth_mixed_jets.py` defines an `AxialSixthMixedJet`
and derives a `RegularFluxFifthMixedJet` analytically. The already-landed
fourth-mixed Eq. (5.2) result remains authoritative for all lower rows.

For the Eq. (5.2) numerator `N = L beta`, with `c = 1/2-h+lambda_n`, the
implementation uses the exact identity

`d_eta^q N_p = 2 eta U_(p,q) + 2 q U_(p,q-1) - (1-eta^2) A_(p,q+1) + 2 eta (q-c) A_(p,q) + q(q-1-2c) A_(p,q-1)`

at `(p,q)=(2,3),(1,4),(0,5)`. Because `L=1-2h eta^2` is quadratic, the
fifth-mixed beta rows close on the landed fourth-mixed rows plus these three
exact numerator derivatives. Production uses neither finite differences nor
sampled division by `X`.

## Compact Lemma 5.2 repair

`background_moment_repair_sixth_mixed_jets.py` lifts the genuine compact
five-bump Lemma 5.2 U repair by one eta tier. On active U-repair support,
`U_XXetaetaetaeta` receives `alpha_j'''' d_X^2 b_j`,
`U_Xetaetaetaetaeta` receives `alpha_j''''' d_X b_j`, and
`U_etaetaetaetaetaeta` receives `alpha_j'''''' b_j`.

The landed fifth-mixed repaired U adapter remains authoritative for every lower
row. Off compact support the unrepaired sixth-mixed jet is returned exactly, so
unneeded sixth moment/patch rows are not requested.

## Hierarchy ownership

`Section5SixthMixedCoefficientJetSource` owns the repaired sixth-mixed U
provider together with the already-landed repaired fifth-mixed phi provider.
`Section5LowerHistorySixthMixedHierarchy` checks exact sixth-to-fifth U
projection coherence and derives `beta_fifth_mixed_jet` only through the
analytic Eq. (5.2) map. No caller-maintained beta derivative table is accepted.

The unrepaired sixth-mixed U provider, unrepaired fifth-mixed phi provider,
normalization `C`, moment/patch eta jets, and genuine order-zero strong profile
remain upstream contracts. This increment does not manufacture the Issue #1
leading profile.

## Verification

`tests/test_background_repaired_history_sixth_mixed.py` uses a real
`Lemma52MomentRepair.from_intervals(...)` five-bump repair. It checks that the
compact repaired sixth-mixed U rows agree with centered eta differentiation of
the landed analytic fifth-mixed repair; that the fifth-mixed Eq. (5.2) beta rows
agree with centered eta differentiation of the landed analytic fourth-mixed
Eq. (5.2) path; that the strong hierarchy owns the repaired sixth U jet and
derives the same fifth-mixed beta as an independent repair adapter; that an
older phi-fifth/U-fifth source fails closed when sixth-mixed U ownership is
requested; and that compact-support short-circuiting does not demand unavailable
sixth eta rows away from repair support. Finite differences occur only as test
oracles. Polynomial jets are analytic test fixtures, not paper coefficient data.

## Truth boundary

Status remains Stage-2 `formal-structure` with `full_reconstruction=false` and
`paper_exact_velocity_available=false`.

This increment supplies the hierarchy-owned fifth-mixed regular-flux input
needed by a future analytic `partial_eta^3(Omega/X)` bridge. It does not itself
establish that third-eta Omega row, the complete `partial_eta^3 actualLowerSource`,
hierarchy-owned `partial_eta^3 f_n`, the third eta jet of `W_n^(0)`, the second
eta jet of `W_n^(1)`, Picard convergence, final coefficient materialization,
uniform `C[j,m]` bounds, recursive cutoff selection, Proposition 5.3 all-order
residual decay, full reconstruction, or paper-exact velocity.

Uniform bounds, recursive cutoff scheduling, and all-order convergence remain
outside this increment and are left to the separate convergence workstream.
