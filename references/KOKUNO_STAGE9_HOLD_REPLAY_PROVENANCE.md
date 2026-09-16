# Kokuno stage-9 hold-scaling replay provenance

## Scope

`KOKUNO5-STAGE9-HOLD-SCALING-REPLAY-002` independently reimplements one exact finite-stage identity displayed in the public Kokuno Navier–Stokes workbench. It checks the endpoint scaling across the outer-schedule hold where `l=-1` for length `4 log(1/h)`.

This increment does **not** copy or execute Kokuno's `proof_sources/stage9/exact_checks.py`. The source repository does not expose a clear repository-level reuse license, so the target checker is independently authored from the public mathematical formulas and records provenance only.

## Pinned public source

- source repository: `KokunoYumeto/yang-mills-interacting-workbench`
- research-state source commit: `e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f`
- public workbench commit inspected: `fab69fdc4ac197159b8e6ae8d73a82bde2b20d55`
- workbench path: `navier-stokes/navier_stokes_workbench.tex`
- workbench Git blob SHA-1: `205a99807302e21a51c5eaf223390c0dfc42bcd0`
- corrected bundle SHA-256 recorded by public metadata: `43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5`
- corrected release DOI: `10.5281/zenodo.22678406`
- component: `NS-stage9`, source pages `100–116`
- archived body path: `proof_sources/stage9/stage9_body.tex`
- archived body SHA-256: `bce636fb3ed7348ed9c73912183f81a75557cc975418ce9e162c4be994f82cb0`
- associated checker path recorded in research-state: `proof_sources/stage9/exact_checks.py`

The bundle/member hash above is recorded provenance. The source ZIP was not extracted in this run, so no fresh byte-level bundle/member re-hash is claimed.

## Independently replayed identity

The public workbench states that on the `l=-1` hold,

- `d/dy log(X E^2) = 2l = -2`,
- `d/dy log(E) = l - 1/2 = -3/2`,
- the hold length is `4 log(1/h)`.

Therefore exact integration gives

- `log((X E^2)_end/(X E^2)_start) = -2 * 4 log(1/h) = 8 log h`, hence the ratio is `h^8`;
- `log(E_end/E_start) = -(3/2) * 4 log(1/h) = 6 log h`, hence the ratio is `h^6`.

The target checker represents the rates, hold-length coefficient, and resulting powers using `fractions.Fraction`. It does not numerically evaluate logarithms. Tolerance is exactly zero.

For the exact rational fixture `h=1/16`, the independently evaluated ratios are

- `(X E^2)_end/(X E^2)_start = 1/16^8 = 1/4294967296`,
- `E_end/E_start = 1/16^6 = 1/16777216`.

A fail-closed regression changes the hold-length coefficient from `4` to `4 + 2^-40`. The resulting exponents are no longer exactly `8` and `6`, and the replay rejects the mutation.

## Actual commands and results

`PYTHONPATH=/tmp/kokuno5_stage9/src python -m pytest -q -W error /tmp/kokuno5_stage9/tests/test_kokuno_stage9_hold_replay.py`

Result: `5 passed in 0.06s`.

`python -m py_compile /tmp/kokuno5_stage9/src/openai_ns_reconstruction/kokuno_stage9_hold_replay.py /tmp/kokuno5_stage9/tests/test_kokuno_stage9_hold_replay.py`

Result: exit `0`.

## Remaining boundary

This verifies only one formula-level stage-9 scaling unit. It does not replay the original Kokuno checker, does not reproduce all recorded stage-9 result groups, and does not establish all-stage realization/convergence, final forcing for the same witness, imported profile existence, or formal endpoint/Comparator closure. A source-recorded PASS is not treated as local evidence.

`paper_exact_velocity_available=false`

`full_reconstruction=false`

Acceptance remains `pending`.
