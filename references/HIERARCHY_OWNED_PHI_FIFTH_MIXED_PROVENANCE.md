# Hierarchy-owned repaired phi fifth-mixed provenance

## Scope

This artifact is one bounded Stage-2 / Issue #2 ownership increment for the
Section 5 positive-order solver.  PR #276 landed the analytic Lemma 5.2 compact
repair for the fifth-mixed angular jet

`phi_n = C E_n / R`.

This layer makes that repaired fifth-mixed `phi_n` provider part of the same
strict lower-history hierarchy that already owns repaired fifth-mixed `U_n`
data.  It does not construct a new profile or infer any coefficient from
samples.

## Exact ownership relation

For each coefficient order with a
`Section5PhiFifthMixedCoefficientJetSource`, one authoritative
`ProfileFifthMixedJet` supplies

- `phi_XXetaetaeta`,
- `phi_Xetaetaetaeta`, and
- `phi_etaetaetaetaeta`,

together with all lower entries.

The inherited fourth-, third-, and second-mixed angular providers must be exact
projections of that fifth-mixed provider.  Every hierarchy query checks this
coherence and fails closed on disagreement.  The existing fifth-mixed `U_n`
provider remains authoritative for the axial hierarchy, so `beta = V/X`
continues to be derived only through the analytic Eq. (5.2) path; no independent
beta derivative table is admitted.

For positive orders, the convenience constructor uses the actual
`Lemma52RepairedPhiFifthMixedJetAdapter` and
`Lemma52RepairedFifthMixedJetAdapter`.  Lemma 5.2 repair is not applied to order
zero; genuine leading strong data remain an Issue #1 dependency and must be
supplied explicitly.

## Verification

`tests/test_background_repaired_history_phi_fifth_mixed.py` checks that

1. the hierarchy-owned repaired fifth-mixed `phi_n` is exactly equal to an
   independently instantiated Lemma 5.2 fifth-mixed adapter;
2. fifth -> fourth -> third -> second angular projections agree exactly with
   the already-landed hierarchy layers;
3. the same hierarchy still owns the repaired fifth-mixed `U_n` jet;
4. an older fourth-mixed-only source fails closed when fifth-mixed `phi_n` is
   requested;
5. an incoherent fifth-to-fourth projection is rejected; and
6. the Lemma 5.2 constructor rejects order zero.

The regression uses the real five-bump `Lemma52MomentRepair.from_intervals(...)`
path.  Polynomial data are test fixtures only.

## Truth boundary

Status remains `formal-structure`.

The unrepaired fifth-mixed `phi_n` provider, fifth-order moment/patch eta jets,
normalization `C`, and genuine order-zero strong profile remain upstream
contracts.  Caller-supplied jets and test polynomials are not paper-exact
coefficient data.

This artifact does **not** establish the third-eta angular preceding-diffusion
row, hierarchy-owned `partial_eta^3 actualLowerSource`, hierarchy-owned
`partial_eta^3 f_n`, `partial_eta^2 W_n^(1)`, Picard convergence, finalized
positive-order coefficient profiles, the paper recursive cutoff-scale schedule,
Proposition 5.3 arbitrary-order residual decay, full reconstruction, or
paper-exact velocity.  `full_reconstruction=false` and
`paper_exact_velocity_available=false` remain mandatory.
