# Hierarchy-owned repaired phi fourth-mixed provenance

## Scope

This increment makes the already-landed compact Lemma 5.2 repaired
`ProfileFourthMixedJet` for `phi_n` part of the authoritative Section 5
lower-history hierarchy.

The new strong coefficient source owns two stronger providers together:

- repaired fourth-mixed `phi_n`, required by the analytic second-eta angular
  preceding-diffusion row;
- repaired fifth-mixed `U_n`, already required by the analytic Eq. (5.2)
  fourth-mixed regular-flux path and the second-eta Eq. (5.6) `Omega/X` row.

The fourth-mixed phi provider is authoritative. Its third- and second-mixed
projections are the inherited phi providers and every hierarchy query checks
that those projections agree exactly. The fifth-mixed U provider remains
authoritative for its lower U projections. No independent beta derivative
table is admitted: beta continues to be derived only through the landed
analytic Eq. (5.2) adapters.

## Paper / implementation map

- Lemma 5.2 and Appendix A: compact five-bump moment repair.
- Eqs. (5.3)-(5.6): recursive strict-lower source whose angular
  preceding-diffusion term needs the stronger phi derivative layer.
- Eq. (5.2): regular-flux beta derivatives remain derived from the owned U
  hierarchy; this increment does not reimplement the `D + lambda_n` shift.
- Eq. (5.6): the existing hierarchy-owned second-eta `Omega/X` path continues
  to consume the inherited fifth-mixed U hierarchy unchanged.
- Eq. (5.7): downstream Picard differentiation will consume the future complete
  hierarchy-owned second forcing derivative.

Implementation artifact:
`src/openai_ns_reconstruction/background_repaired_history_phi_fourth_mixed.py`.

## Verification

Regression uses the actual `Lemma52MomentRepair.from_intervals(...)` compact
five-bump geometry and analytic moment/patch eta jets.

1. The hierarchy-owned repaired fourth-mixed phi jet is compared exactly with
   an independently instantiated `Lemma52RepairedPhiFourthMixedJetAdapter`.
2. Its third- and second-mixed projections are required to equal the existing
   hierarchy-owned lower phi layers exactly.
3. The same hierarchy is checked to preserve the already-landed repaired
   fifth-mixed U ownership rather than creating a parallel strong-data table.
4. The owned fourth-mixed phi jet is passed through the already-landed analytic
   second-eta preceding-diffusion operator and compared with an independently
   instantiated repaired phi provider.
5. A previous fifth-mixed-U/third-mixed-phi source fails closed when the new
   fourth-mixed phi path is requested.
6. A deliberately incoherent fourth-to-third phi projection is rejected.
7. The Lemma 5.2 factory remains positive-order only; order zero must be
   supplied explicitly once genuine leading strong data are available.

No production finite difference, sampled `Omega/X`, generic cutoff, or fitted
coefficient family is introduced.

## Truth boundary

Status remains **Stage 2 / `formal-structure`** and
`paper_exact_velocity_available=false`.

The unrepaired/base fourth-mixed phi jet, unrepaired/base fifth-mixed U jet,
normalization `C`, moment/patch eta jets, and especially the genuine order-zero
strong profile remain upstream inputs. Issue #1 still owns the real leading
strong data. Test polynomials and explicit caller providers are validation
fixtures, not paper-exact coefficients.

This increment establishes only hierarchy ownership and coherence of the
strong phi fourth-mixed layer. It does **not** yet claim a complete
hierarchy-owned `partial_eta^2 actualLowerSource`, `partial_eta^2 f_n`,
`partial_eta W_n^(1)`, Picard convergence, final positive-order coefficient
materialization, the recursive cutoff-scale schedule, Proposition 5.3 all-jets
residual decay, or a paper-exact reconstructed velocity.
