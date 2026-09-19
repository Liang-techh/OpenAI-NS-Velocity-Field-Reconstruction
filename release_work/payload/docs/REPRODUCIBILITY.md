# Reproducibility

## Python

From the repository root:

```bash
python -m pip install -r requirements.txt
python ns_candidate.py verify
python ns_candidate.py evaluate --candidate ST054-Q2 --point 0.1 0 0.1 --time 0.5
python -m pytest -q -W error tests
python ns_candidate.py validate --candidate ST054-Q2 --seed 9175491 --out outputs/Q2_recheck.json
```

The last command is expected to exit 1 because the momentum gates fail. It refuses to overwrite an existing file. Repeat with ST054-M3 or seed 9175492. Full validation uses 4096 Cartesian points, six times, three spatial steps, separate temporal refinement and three energy-quadrature orders. No training or network access is needed after dependencies are installed.

The repository-local API is `from ns_candidate import load`; `family, raw = load('ST054-Q2')`, then `family.fields(raw, points, time)` or `family.analytic_residual(raw, points, time)`. Use batches of 256-1024 points for large derivative evaluations. This preserves the original evaluator rather than introducing another optimization model.

## Identity

The native-tested MAT contains frozen raw coefficient vectors and Python reference values. Release JSON candidates are serialized from those vectors with new English metadata. Their JSON hashes intentionally differ from original research JSONs. The registry records both identities. Reference u/p/full-residual values and original validation norms must match documented tolerances.

The runtime basis/evaluator files and independent validator are copied from a pinned source commit. The MAT is copied unchanged from its pinned native-tested publication, with fixed SHA256. manifest.json binds packaged evidence, code and data. Temporary outputs are outside that inventory.

## MATLAB

```matlab
addpath('visualization/matlab');
ns_explorer;
ns_core_explorer;
ns_release_test('outputs/matlab');
```

Normal viewing needs MATLAB and the included MAT, not Python. The source app was tested with R2024b. This release's own MATLAB workflow checks its evaluator, sliders, playback and core tools; only completed results establish what ran. A source receipt is historical, not a claim that a new release ran those tests.

Legacy tests for removed reconstruction modules remain in backup history. No full legacy suite or Lean build is claimed. Software/reference tolerances around 1e-8 must not be mistaken for the actual nonzero NS residual around 0.04.
