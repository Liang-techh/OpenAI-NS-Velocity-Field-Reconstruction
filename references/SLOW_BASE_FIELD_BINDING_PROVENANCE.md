# Slow-box / theorem physical-base binding provenance

Status: **formal-structure only**.

Base main: `1d812ff01d4575551f804c3599d6efa8985c7d2f`.

This increment closes one Section 6->7 identity seam without inventing the
still-missing paper-exact field values. The executable slow-label path already
accepts a `TangentialBaseJetProvider`, while the current main separately
contains a theorem-facing `LargeBandPhysicalBaseBinding` which identifies one
admitted slow box with the pinned
`PrimaryGeometryAssembly.frequency_eq_physical` /
`PrimaryGeometryAssembly.axial_eq_physical` chain ending at
`FinalSlowBase.velocity`.

`slow_base_field_binding.py` now makes those two identities composable. A
provider registry entry must name the exact sign-free slow box, source revision,
and selected `Prepared` instance used by the theorem physical-base binding
before it can be described as theorem-backed.

## Inputs inspected on the pinned base

- `README.md` and `AGENTS.md`, including the explicit prohibition on promoting
  diagnostic or theorem-shaped plumbing to a completed reconstruction.
- `references/provenance_manifest.json`, whose official Lean pin remains
  `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
- `src/openai_ns_reconstruction/charts.py`: binary64 Eqs. (6.1)-(6.6)
  geometry, with runtime `DyadicChart` restricted to `ell<=1000`.
- `src/openai_ns_reconstruction/slow_labels.py`: executable labels and the
  caller-supplied tangential-base-jet provider path, likewise binary64-band
  bounded.
- `src/openai_ns_reconstruction/phase.py`: local Eqs. (7.1)-(7.8) only.
- `src/openai_ns_reconstruction/phase_large_band_scale.py`: theorem-facing
  large-band scalar phase gate that intentionally does not form binary64 `Q`.
- `src/openai_ns_reconstruction/phase_large_band_prepared_source.py`: pinned
  `PrimaryGeometryAssembly.Prepared` source identity.
- `src/openai_ns_reconstruction/phase_large_band_physical_base.py`: pinned
  theorem-facing bridge from that Prepared source to components of
  `FinalSlowBase.velocity`.

Relative to the previous tested head, current main adds merged PR #283, a
bounded Section-5 hierarchy-owned third-eta axial preceding-diffusion increment.
`aa3f0ee... -> 1d812ff...` changes only its Section-5 source/test/provenance
files, so it does not overlap this four-file Sections 6-7 increment. Open PRs
#284 and #286 are Stage-1 work and likewise do not overlap. No explicit current
`Agent 8`-named claim is visible in the accessible Issue/open-PR records, but
Sections 8-9 already have substantial parallel mean-correction/residual-iteration
formal layers; this increment therefore treats them as externally owned and
makes no edits there.

## What this increment verifies

The `SlowBoxBaseFieldBindingWitness` remains sign-free and still records an
opaque normalized-field/provider identity. It now additionally records the
`source_id` and `prepared_instance_id` needed to compare against the landed
theorem-facing base source.

`TheoremBackedSlowBoxBaseFieldBinding` fails closed unless all of the following
are literally identical on both sides:

1. the unsigned slow-box key `(ell,a)`;
2. the `(source_id, source_revision)` pair;
3. the `(source_id, source_revision, prepared_instance_id)` Prepared key; and
4. the already-landed physical-base theorem binding remains fully certified.

A successful composition therefore says only that the provider registry points
at the same formal slow base already identified with the pinned
`FinalSlowBase.velocity` chain. It does not evaluate that field or certify the
provider's returned numerical jet.

## Large-band / binary64 boundary

The theorem-facing Section-7 scale adapters are eventual large-band statements
and may name `ell>1000`; the current executable `DyadicChart` / `SlowLabel` API
is intentionally restricted to `ell<=1000`. This increment therefore allows
the *metadata identity* to carry any positive integer band while exposing
`binary64_runtime_box_supported`. It does not coerce a deep asymptotic band
into the binary64 runtime and does not sample a surrogate field to bridge that
gap.

The ordinary `freeze_bound_label_phase_data` path still requires an actual
`ActiveSlowRepresentative`, so its executable band boundary is unchanged.

## Verification

`tests/test_slow_base_field_binding.py` keeps the executable sign-pair/provider
regressions and adds an independent theorem-facing fixture assembled from the
landed `LargeBandPhaseScaleCertificate`, `LargeBandLocalBaseAdmission`,
`LargeBandBaseSourceBinding`, `LargeBandPreparedSourceBinding`, and
`LargeBandPhysicalBaseBinding` types. It verifies exact identity composition,
cross-wiring rejection, and the explicit deep-band/binary64 boundary. The
fixture uses synthetic theorem metadata and is not manuscript field data.

The identical source/test payload passed complete repository Actions run #1997
on the immediately preceding `aa3f0ee...` base, including both full Python/NumPy
matrices, wheel/outside-checkout CLI/audit, paper-exact refusal, and the four
slices. Because the base has advanced to `1d812ff...`, this rebased head still
requires its own current-base workflow before merge.

## Truth boundary

This module does **not** replay `canonicalPrepared`, materialize
`FinalSlowBase.velocity`, derive a real normalized tangential jet, certify
LocalBase C1/C2 values, prove the all-active-family applications of
Eqs. (7.9)-(7.11), solve the primary-amplitude ODE on paper data, or construct a
supported divergence-free paper wave.

Accordingly:

- `actual_physical_base_values_materialized = false`
- `base_jet_values_machine_certified = false`
- `uniform_phase_estimates_proved = false`
- `paper_exact_velocity_available = false`

No surrogate background, profile, amplitude, or wave is introduced.
