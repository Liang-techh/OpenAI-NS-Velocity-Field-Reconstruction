# Navier-Stokes Candidate Lab

Reproducible finite-window velocity fields, independent full-residual validation, and interactive MATLAB visualization.

**Two complete frozen research candidates: ST054-Q2 (lower L2) and ST054-M3 (lower peak). The original full NS 1e-3 target is not yet met.** This is an independent project, not an OpenAI repository or a reproduction of its exact field.

![Native MATLAB field explorer](evidence/source_matlab/matlab_explorer.png)

*Actual native MATLAB screenshot from the pinned source publication. Full-field and dense central-core viewers are included; visual similarity does not certify the PDE.*

## Start with the visualization

Download or clone this repository, open MATLAB in its root folder, and run:

```matlab
addpath('visualization/matlab');
ns_explorer;
```

Drag the bottom physical-time slider to evaluate the frozen field continuously over **0.25 to 0.75**. Rotate the 3-D view, move slices, change seed positions, compare pressure force, vorticity and full residual, or switch candidates. Normal viewing needs only MATLAB and the bundled data, not Python, training or a network download.

For high-resolution local sampling rather than a cropped full-domain grid:

```matlab
ns_core_explorer;
T = ns_core_timeseries([], 'ST054-Q2', 'Times', linspace(0.25,0.75,13));
```

[Visualization guide](visualization/matlab/README.md) | [Scientific status](docs/SCIENTIFIC_STATUS.md) | [Reproducibility](docs/REPRODUCIBILITY.md)

## Frozen independent results

Each row summarizes 4096 Cartesian points, six fixed times and the original refinement ladders. L2 is spatial volume L2, not RMS. Compare candidates on the same seed.

| Candidate | Seed | Sampled maximum | Spatial volume L2 |
|---|---:|---:|---:|
| ST054-Q2 | 9175491 | 0.03730219 | 0.04761902 |
| ST054-Q2 | 9175492 | 0.03845539 | 0.04770423 |
| ST054-M3 | 9175491 | 0.03640765 | 0.04832320 |
| ST054-M3 | 9175492 | 0.03712486 | 0.04789221 |

Both momentum gates remain **FAILED**. Other recorded numerical gates pass on those samples; there is no continuous-domain error bound or blow-up proof. The release replays are in [evidence/replay](evidence/replay), with original source summaries in [evidence/source_results.json](evidence/source_results.json). The lower-L2 default is a publication choice among these two complete candidates, not a claim of global optimality.

## Evaluate and verify in Python

```bash
python -m pip install -r requirements.txt
python ns_candidate.py verify
python ns_candidate.py evaluate --candidate ST054-Q2 --point 0.1 0 0.1 --time 0.5
python -m pytest -q -W error tests
python ns_candidate.py validate --candidate ST054-Q2 --seed 9175491 --out outputs/Q2_recheck.json
```

The final command intentionally exits **1** while the scientific gates fail. Outputs are never silently overwritten. Software CI success means reproducibility, not PDE acceptance.

## Repository map

| Directory | Purpose |
|---|---|
| `candidates/` | Complete Q2/M3 coefficients, immutable identities and selection scope |
| `runtime/` | Pinned finite-basis evaluator and original independent Cartesian validator |
| `visualization/matlab/` | Full-field explorer, dense core explorer, time-series diagnostics and model MAT |
| `evidence/` | Source results, four full replays and numerical reference parity |
| `tests/` | Focused release regressions, separate from scientific acceptance |
| `docs/` | Method boundaries, reproduction and migration record |

## Provenance and previous work

This is a clean publication of completed results from [the research workspace](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1), not a wholesale merge of unfinished experiments. Source commits and SHA256 identities are pinned in [provenance.json](provenance.json). Candidate JSON metadata is newly serialized from the immutable native-tested coefficient bundle; original research and release file hashes are explicitly distinguished.

The previous repository tree is preserved in `archive/pre-st054-replacement-2026-09-19` and Git history. Old issues, PRs and other branches are not deleted. ST055 intermediate comments without a complete candidate package are not promoted. See [the migration record](docs/MIGRATION.md).

`pde_validated=false` | `source_correspondence_verified=false` | `paper_exact=false` | `blowup_proved=false`
