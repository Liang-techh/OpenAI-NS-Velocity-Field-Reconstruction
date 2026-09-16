# Kokuno base-heat replay provenance

## Scope

`KOKUNO5-BASE-HEAT-ODE-REPLAY-001` independently reimplements one exact algebraic sub-check from the public Navier–Stokes base/heat reconstruction. It checks the polynomial reduction used to derive the heat-factor ODE and the matching total-derivative factor, plus the adjacent exact coefficient identity `2 A^2 - 1/2 = 2 h(1+h)` for `A = 1/2 + h`.

This is **not** a copy or execution of the Kokuno checker. It is a clean reimplementation from the mathematical formulas displayed in the public workbench.

## Pinned public source identities

- source repository: `KokunoYumeto/yang-mills-interacting-workbench`
- published research-state source commit: `e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f`
- public workbench commit inspected for the displayed formulas: `fab69fdc4ac197159b8e6ae8d73a82bde2b20d55`
- `navier-stokes/navier_stokes_workbench.tex` Git blob SHA-1: `205a99807302e21a51c5eaf223390c0dfc42bcd0`
- corrected source bundle SHA-256 recorded by the public metadata: `43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5`
- corrected release DOI: `10.5281/zenodo.22678406`
- base-heat derivation member SHA-256: `98d101542a7d8d5c1475764ec99178a663fe363e416d572eeafa048e724a193d`
- recorded base-heat checker path: `proof_sources/base_heat/verify_identities.py`
- recorded checker member SHA-256: `da6fe495af223772824b2f658c83c4aec08d7aa41882818589eaf14c6a3425dc`
- recorded result path: `proof_sources/base_heat/identity_check_results.json`
- recorded result member SHA-256: `5b74d2d2598841ce2512790dd70efe35640b69a4ea5a8516a4dbbd482f0bd843`

The bundle/checker/result SHA-256 values above were cross-checked against `navier_stokes_checks.json` metadata. The binary source ZIP was not available through the GitHub text connector in this run, so these member hashes were **not** independently recomputed from extracted bytes.

## License / copyright boundary

No repository-level `LICENSE`, `LICENSE.md`, `COPYING`, or `NOTICE` grant was found for the Kokuno repository, and the corrected Zenodo record does not expose a usable license value in the checked metadata. Therefore this increment copies no Kokuno Python checker code, proof prose, or result JSON. The implementation below is independently authored from the public mathematical identities and retains URL/commit/hash provenance.

## Independently replayed identity

The public workbench defines the heat factor

`H(Z) = Gamma(1+h)^(-1) integral_0^infinity exp(-v) v^h (1+Zv)^(-h) dv`

and, after putting the ODE terms over one common integrand, displays the exact polynomial reduction

`(1+h) Z^2 v^2 - (1 + 2(1+h)Z) v(1+Zv) + (1+h)(1+Zv)^2`

`= (1+h) - v - Z v^2`.

It further identifies the right-hand side as the factor obtained from

`d/dv [exp(-v) v^(1+h) (1+Zv)^(-1-h)]`

after removing the same common factor. The boundary values at `v=0` and `v=infinity` are stated as zero, yielding the heat-factor ODE

`Z^2 H'' + (1 + 2(1+h)Z) H' + h(1+h) H = 0`.

The target replay expands the two polynomial expressions independently over the exact rational polynomial ring `Q[h,Z,v]`. No floating-point evaluation, sampling, quadrature, or numerical tolerance is used.

## Actual replay result

Exact reduced coefficient map, with monomials ordered as `(h-power, Z-power, v-power)`:

- `(0,0,0) -> 1`
- `(1,0,0) -> 1`
- `(0,0,1) -> -1`
- `(0,1,2) -> -1`

Both the ODE common-factor expression and the total-derivative expression reduce to exactly that map. The radial coefficient identity also reduces exactly on both sides. Tolerance is exactly `0`.

Failure-boundary regression perturbs the constant coefficient by `2^-40`; exact coefficient comparison rejects the perturbed polynomial. Thus there is no near-zero acceptance threshold in this checker.

Commands actually executed against the authored files before push:

`PYTHONPATH=/tmp/kokuno5_replay/src python -m pytest -q -W error /tmp/kokuno5_replay/tests/test_kokuno_base_heat_replay.py`

Result: `3 passed in 0.06s`.

`python -m py_compile /tmp/kokuno5_replay/src/openai_ns_reconstruction/kokuno_base_heat_replay.py /tmp/kokuno5_replay/tests/test_kokuno_base_heat_replay.py`

Result: exit `0`.

## Remaining unverified scope

- The original `proof_sources/base_heat/verify_identities.py` bytes were not extracted or executed.
- The original `identity_check_results.json` was not treated as local evidence of PASS; only its published member hash was recorded.
- The rest of the base-heat checker/result groups are not reproduced by this increment.
- None of the source record's nine replay programs is claimed fully reproduced by this one sub-check.
- The other recorded 144 named mathematical result groups and 12 structural probes remain outside this increment except for this independently named target-side result group.
- This does not validate imported profile existence, all-stage realization/convergence, final forcing, formal endpoint/Comparator closure, or the complete theorem.

`paper_exact_velocity_available=false`

`full_reconstruction=false`
