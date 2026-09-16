# Kokuno Navier–Stokes source import map

Status: provenance-only migration map; no source proof body or source code copied.

Target baseline: `codex/stage1-complete-slow2` at `77ed17c1e81b7abce99fa4c7ffc5b3d96ba389d0`.

Source repository: `https://github.com/KokunoYumeto/yang-mills-interacting-workbench`, `navier-stokes/` public record. The source record identifies `e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f` as its `source_commit`; that Git commit exists and is titled `Organize Yang-Mills, S6, and Navier-Stokes research`. Corrected archived edition: DOI `10.5281/zenodo.22678406` (9 September 2026).

This document inventories facts that may safely be carried into this repository and maps the thirteen published source components onto existing reconstruction tasks. It is **not** an independent mathematical validation, a paper-exact certificate, or permission to copy the source bodies.

## Truth and copyright boundary

The Kokuno repository is publicly readable, but this audit found no root `LICENSE`, `LICENSE.md`, `COPYING`, or `NOTICE` file. The corrected Zenodo record exposes a `Rights / License` heading but no license value. Public availability alone is therefore not treated here as a copyright license.

Operational import policy for this repository:

- **DIRECT-METADATA**: factual identifiers may be recorded: repository/DOI/URL, commit IDs, archive/member paths, file sizes, hashes, check counts, and independently authored descriptions of relationships.
- **REIMPLEMENT**: source mathematical prose, LaTeX bodies, Python checkers, and other copyrightable source content are not copied verbatim or in bulk from the Kokuno bundle while its reuse license remains unclear. Required functionality should be independently implemented from the cited primary/formal sources and checked against factual hashes/provenance.
- **REFERENCE-ONLY**: source review/audit conclusions may guide what to inspect, but are not imported as proof evidence.
- **UNVERIFIED**: any claim not independently replayed here remains explicitly unverified.

The separately cited formal repository `openai/NavierStokesAndEuler` at commit `8937a8f4cbc7abaab5e9e97d1cc7f5d2319d9538` contains an Apache-2.0 `LICENSE`. That licensing fact applies to that formal repository, not automatically to Kokuno-authored reconstruction text, bundle members, or the primary manuscript. If formal-source code is later reused directly, Apache-2.0 attribution/redistribution requirements must be followed separately.

No material copied by this increment relies on an assumption that the primary manuscript itself grants a redistribution license.

## Source identities and hash inventory

The following values are recorded by the public source metadata. `navier-stokes/research-state.json`, `navier-stokes/navier_stokes_checks.json`, and `releases/2026-09-09/HASHES.md` agree on the corrected source ZIP SHA-256. The source checks also record `member_count=120`, `embedded_text_files=118`, `all_member_hashes_verified=true`, `all_component_paths_present=true`, and a successful fresh extraction/replay/rebuild. Those are **source-recorded check results**; this run did not independently recompute the ZIP or all 120 member hashes.

| Item | Recorded identity | Import status |
| --- | --- | --- |
| source record commit | `e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f` | DIRECT-METADATA; commit existence checked |
| primary manuscript capture used by component/profile/stage derivations | 165 pages; SHA-256 `8c8a94ad9ac824c8b605b9827cadf7beaca48bd10b380de3cfc872a2c37afa81` | DIRECT-METADATA; bytes not rehashed here |
| later manuscript capture used additionally by summation/radial/cycle audits | 166 pages; SHA-256 `0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f` | DIRECT-METADATA; bytes not rehashed here |
| formal source | `https://github.com/openai/NavierStokesAndEuler/tree/8937a8f4cbc7abaab5e9e97d1cc7f5d2319d9538` | DIRECT-METADATA; Apache-2.0 repository license observed |
| corrected source bundle | `navier_stokes_source_bundle.zip`; SHA-256 `43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5`; 120 members | DIRECT-METADATA only until local byte rehash is possible |
| corrected reader PDF | SHA-256 `242fe7f83feab43b91915595230b5e93e9f344eb654e1edcb65e9669a9cd5878` | DIRECT-METADATA |
| source check JSON identity | SHA-256 `386d3fef8edcf2665ecd5daa2944269f414ced12c1afe89d9d7fc50973fbd76a` | DIRECT-METADATA |
| corrected Zenodo edition | DOI `10.5281/zenodo.22678406` | DIRECT-METADATA; no displayed license value observed |

The public release map additionally lists the corrected Navier–Stokes ZIP under the same SHA-256 `43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5`. The Zenodo file listing names `released-ns-reader-corrected-2026-09-09.zip` and the 4.4 MB `released_ns_primary_corrected_20260909.json`; file naming differences must not be treated as byte identity without the release mapping/hash evidence.

### Component SHA-256 inventory

These are the thirteen public component identities recorded by `research-state.json`. The source package reports that all member hashes were verified; this target audit records the identities but does not elevate that source-side check to independent validation.

| Component ID | Bundle member | Recorded SHA-256 |
| --- | --- | --- |
| `NS-profiles` | `proof_sources/profiles/profiles_body.tex` | `63327f4a6d339d39de096230af1810cb83e7c427ba0b0cc78c576b2750437a95` |
| `NS-base_heat` | `proof_sources/base_heat/derivation.tex` | `98d101542a7d8d5c1475764ec99178a663fe363e416d572eeafa048e724a193d` |
| `NS-oscillations` | `proof_sources/oscillations/oscillations_body.tex` | `3fd5c61de39b6f65a581ead5b29c741c2e0f46fb30f62e582dc24d93bbe83615` |
| `NS-mean_corrections` | `proof_sources/mean_corrections/mean_corrections_body.tex` | `c494cebcaf476a01c5aae17d38f84dfc6b0092040c03b9aecbad1f5c3831ea29` |
| `NS-stage9` | `proof_sources/stage9/stage9_body.tex` | `bce636fb3ed7348ed9c73912183f81a75557cc975418ce9e162c4be994f82cb0` |
| `NS-terminal_endpoint` | `proof_sources/terminal_endpoint/endpoint_derivation.tex` | `ab3a4365baabf9203aeeb5ea2c74c4cbdcfe897dd9198aa680e19d7178e16358` |
| `NS-profile_assembly` | `proof_sources/profile_assembly_audit/profile_assembly_body.tex` | `32ed054afdb0aca39ad66aa03a30ef9a3d3e004b60d7d9ba1b972435703c3c28` |
| `NS-stage_inputs` | `proof_sources/stage_inputs_audit/stage_inputs_body.tex` | `61aa4684591b964321b00304d58ca7e8e2bbf538cfc930bc8a6be219755bd7e3` |
| `NS-summation_realization` | `proof_sources/primary_continuation/summation_realization_body.tex` | `13d370917df47e46e2d75e0cd3a5d4ab5f2dfd1fd1f31e4eef3fa4413b6a3fa5` |
| `NS-terminal_trace` | `proof_sources/primary_continuation/terminal_trace_body.tex` | `f0eb87504ca5246eec37b09cfb17807053eee37e5282bc6ea4badefb3c822eb5` |
| `NS-radial_stress_reconstruction_audit` | `proof_sources/lanes/web_bundle_audit/tex_handoff/portable_upstream/radial_stress_reconstruction_audit.tex` | `dd6c119670316f1804226723caaefe062a580e40d2101efab86d823b5858ed8a` |
| `NS-stage_cycle_gain_audit` | `proof_sources/lanes/web_bundle_audit/tex_handoff/portable_upstream/stage_cycle_gain_audit.tex` | `497cb07b95a651bd5ed76f74c004722dc2172c1ef0575ddfea43c95cd56b361e` |
| `NS-global_closure_body` | `proof_sources/global_analytic_closure_review/global_closure_body.tex` | `b212819c6788fc74086e501292567ba7630cd67f7a308f54c2102732c40e0954` |

## SOURCE_IMPORT_MAP

`Body policy` below is deliberately conservative. Every row permits metadata import, but no Kokuno mathematical body is copied while the Kokuno licensing status is unresolved.

| Source component | Target tasks / existing target surfaces | Existing work to avoid duplicating | Body policy | Current verification status |
| --- | --- | --- | --- | --- |
| `NS-profiles` | NS016 actual leading-profile assembly; NS017 Stage-1 support/moment/matching/cone; `NaturalProfileAssembly` / `LeadingProfile` | open PR #364 exact-average callback; open PR #325 profile-assembly contract | REIMPLEMENT from primary/formal sources; source body REFERENCE-ONLY | hash/path inventoried; mathematical body not independently replayed by this import |
| `NS-base_heat` | NS016–NS017 regular-inner-core / heat-exterior matching and Stage-1 theorem prerequisites | no dedicated checkpoint task PR identified; do not create a parallel profile stack | REIMPLEMENT; factual formulas/section pointers may be cited | hash/path inventoried; source-side identity checks only |
| `NS-oscillations` | NS028 phase/wave-vector/polarization; NS029 stress-cone decomposition and oscillatory-wave producer | NS027 precursor overlaps open PR #361 Prepared/LocalBase; later oscillation implementation must consume actual upstream objects | REIMPLEMENT; source derivation REFERENCE-ONLY | hash/path inventoried; no actual target oscillatory producer certified here |
| `NS-mean_corrections` | NS030 mean correction and five-equation defect solve | existing Section-8 algebra/support modules should be extended, not replaced | REIMPLEMENT | hash/path inventoried; actual mean-correction field remains unverified here |
| `NS-stage9` | NS031 actual one-step residual improvement; NS032 infinite-stage common scale / Eq. 9.21 field | open PRs #328, #334, #338, #343, #350 cover finite/theorem/provenance seams; do not re-create them | REIMPLEMENT actual stage producer/sum; use source audit as REFERENCE-ONLY | source reports finite checks, but actual target infinite stage field remains unverified |
| `NS-terminal_endpoint` | NS033 final space/time localization; NS034 forcing and all-order endpoint smoothness; NS035 final invariants | open theorem-side PRs #329, #335, #355, #360, #363 and energy PR #346 overlap claims/interfaces | REIMPLEMENT actual runtime chain; audit conclusions REFERENCE-ONLY | endpoint/formal replay explicitly unfinished in source record |
| `NS-profile_assembly` | NS016–NS017 same-source profile assembly and theorem prerequisites | open #325 and #364 already cover narrow contracts; do not import a second assembly model | REIMPLEMENT using existing target objects | hash/path inventoried; source assembly audit is not target runtime proof |
| `NS-stage_inputs` | NS027 actual background→Prepared/LocalBase; NS028 phase inputs; NS029 stress/wave inputs | #361 is a formal Prepared/LocalBase binding precursor; consume rather than duplicate if integrated | REIMPLEMENT producer from target hierarchy; source input ledger REFERENCE-ONLY | hash/path inventoried; no source body copied |
| `NS-summation_realization` | NS025 hierarchy-owned all-order bounds; NS026 infinite diagonal cutoff/background; NS032 correction-sequence sum | legacy open finite/convergence lane #290/#303/#312/#315/#318/#324/#327/#333/#337/#342/#345/#349/#354/#358/#362 overlaps theorem/finite-prefix work | REIMPLEMENT actual total provider/infinite realization; source body REFERENCE-ONLY | 166-page source capture also relevant; no independent infinite-sum validation from this import |
| `NS-terminal_trace` | NS033 localization and NS034 endpoint force/smoothness | endpoint/theorem PRs #329/#335/#355/#360/#363 overlap interfaces | REIMPLEMENT actual trace/force chain | hash/path inventoried; source says fresh Lean/endpoint axiom work unfinished |
| `NS-radial_stress_reconstruction_audit` | NS029 stress reconstruction; NS030 mean-defect coupling | no new parallel stress solver should be created without checking existing Section 7–8 modules | REIMPLEMENT; audit text REFERENCE-ONLY | 166-page source capture relevant; hash inventoried, math not independently validated |
| `NS-stage_cycle_gain_audit` | NS031 residual-improvement cycle; NS032 common scale/summed field | #328/#334/#338 already encode finite residual/gain arithmetic seams | REIMPLEMENT actual stage cycle; retain source hash/provenance as comparison | 166-page source capture relevant; source-side finite checks are not all-order target proof |
| `NS-global_closure_body` | NS033–NS036 final localization, forcing, invariants, reproducible audit | open #335/#346/#355/#360/#363 contain theorem-shaped pieces; NS036 remains the eventual independent completion audit | REFERENCE-ONLY until all target dependencies are actual; reimplement/check every required closure in target | source itself states complete independent analytic/formal validation unfinished |

## What can be imported now

Safe immediate imports are limited to provenance facts and machine-readable inventory: source URLs, DOI, source/formal commits, manuscript and archive hashes, member paths/hashes, component IDs, recorded check counts, and a mapping from those identities to target tasks. These facts may be used to fail closed when a future artifact is expected to match a named source edition.

Direct copying of Kokuno proof bodies/checker code is **not authorized by this map**. If a later audit finds an explicit license covering particular files, those rows can be upgraded file-by-file with the exact license source and obligations. Until then, implementation work should be independently authored against the primary manuscript and/or the separately Apache-2.0 formal source, with provenance recorded.

## Verification limits and source-side checks

The corrected source record reports nine replay programs, 144 named mathematical result groups, 12 structural probes, 120 archive members, successful member-hash/path verification, fresh ZIP extraction replay, and reader rebuild. It also explicitly says complete independent analytic and Lean validation is unfinished. This map preserves that boundary.

This run did **not**:

- recompute SHA-256 over the downloaded 4.3 MB ZIP or all 120 members (binary download was unavailable in the current tool path);
- independently rebuild the 208-page reader;
- rerun the nine source replay programs;
- rerun Lean at the recorded formal commit;
- perform a full semantic parse of the ~4.3 MB `navier_stokes_primary_manifest.json` through the GitHub text connector (the file identity was reachable, but the oversized body was not returned usefully);
- run target `pytest`, `ns-reconstruct demo`, `ns-reconstruct audit`, or Lean, because this increment changes documentation/provenance only.

Consequently, hash agreement in this document means **public-metadata cross-check**, not local byte verification. A future byte-level import increment should download the DOI artifact, recompute the archive hash and every member SHA-256, compare the exact member set, and record the commands/results before any file is admitted as byte-verified.

## Actual checks performed for this increment

- target branch lookup `GET /repos/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction/branches/codex/stage1-complete-slow2` → `77ed17c1e81b7abce99fa4c7ffc5b3d96ba389d0`;
- target code search for `SOURCE_IMPORT_MAP` and Kokuno source-import references → no existing target map found;
- target open-PR scan → no open PR dedicated to Kokuno source inventory/license/hash migration; existing mathematical and reconciliation PRs were treated as occupied work;
- source reads: `navier-stokes/README.md`, `RESEARCH_STATE.md`, `research-state.json`, `navier_stokes_checks.json`, directory metadata for `navier_stokes_primary_manifest.json` / source ZIP, plus `releases/2026-09-09/HASHES.md`;
- source commit lookup for `e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f` → commit exists;
- Kokuno license probes `LICENSE`, `LICENSE.md`, `COPYING`, `NOTICE` → all returned 404; repository search found no MIT/Apache/GPL/copyright license declaration;
- formal-source `LICENSE` at `openai/NavierStokesAndEuler@8937a8f4cbc7abaab5e9e97d1cc7f5d2319d9538` → Apache License 2.0;
- corrected Zenodo record → DOI and file inventory present; `Rights / License` rendered without a license value;
- `pytest`: not run (documentation/provenance-only change);
- `ns-reconstruct demo`: not run;
- `ns-reconstruct audit --require-paper-exact`: not run;
- Lean build: not run.

## Mandatory target truth status

This provenance map changes no mathematical completion flag. In particular:

- `paper_exact_velocity_available=false`
- `full_reconstruction=false`
- complete independent validation remains unavailable
- all affected NS task acceptance states remain whatever the coordinator records; this document does not mark any NS001–NS036 task accepted or complete.
