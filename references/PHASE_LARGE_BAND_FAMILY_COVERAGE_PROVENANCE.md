# Section 7 large-band finite active-family coverage provenance

Status: **formal-structure**.

This increment closes one bookkeeping gap after the per-family large-band
`phase_estimates`, `base_estimates -> coordinate_errors`, and `damping_error`
admission bridges. Those gates certify theorem metadata one signed family at a
time; by themselves they do not prevent a caller from proving only a sparse
subset and then describing the resulting bounds as uniform over an active
family.

`phase_large_band_family_coverage.py` therefore accepts an externally supplied,
theorem-provenanced declaration of the **exact finite active signed family set
for one band and one base-source revision**. It then requires exact set
coverage, with no missing, extra, or duplicate labels, from all three existing
theorem paths. For each signed label the coordinate and damping paths must
reuse the same upstream family admission as the phase path. The declared set
is also checked mechanically for the Section 6 `sigma = +/-1` sign-pair
closure.

The formal theorem identities reused by this coverage layer are pinned to:

- repository: `openai/NavierStokesAndEuler`
- commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- module: `NavierStokes/BasePhaseGeometry.lean`
- `BasePhaseGeometry.FamilyData.phase_estimates`
- `BasePhaseGeometry.FamilyData.base_estimates`
- `BasePhaseGeometry.FamilyData.coordinate_errors`
- `BasePhaseGeometry.FamilyData.damping_error`

The upstream formal source states `phase_estimates` for an arbitrary family
index `i`, carrier point `q`, and slot variable `v`, under its large-band and
membership hypotheses. This Python increment does **not** enumerate those
indices from the paper. Instead, exactness of the finite scope is a separate
analytic/formal theorem witness and sampled, fitted, or numeric-scan evidence
is rejected.

A successful coverage object may report maxima of the already-admitted theorem
envelopes over that declared finite scope. Those maxima are not measured field
errors and they are not a global all-band result. In particular, this
increment does not provide the genuine Section 6 active carrier/slot family,
the Proposition 5.5 paper-exact base/background provider, actual phase/base/
frame/damping fields, or a proof that the supplied tuple is the paper's true
active set.

Truth boundary remains:

- `exact_active_family_scope_machine_verified = false`
- `theorem_witnesses_machine_verified = false`
- `global_all_band_uniformity_verified = false`
- `uniform_eq_7_9_to_7_11_verified = false`
- `paper_exact_velocity_available = false`

No surrogate wave, fitted background, or sampled uniformity claim is introduced.
