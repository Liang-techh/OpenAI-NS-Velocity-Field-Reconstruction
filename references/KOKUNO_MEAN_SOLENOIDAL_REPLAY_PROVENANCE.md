# Kokuno NS mean-correction solenoidal replay provenance

Status: `formula-level-replay`; acceptance remains `pending`.

## Source boundary

This increment independently reimplements one short public identity from the
Kokuno Navier--Stokes workbench Section 8 / `NS-mean_corrections` component:

- public workbench repository commit inspected: `fab69fdc4ac197159b8e6ae8d73a82bde2b20d55`
- research-state `source_commit`: `e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f`
- source pages recorded by research-state: `87–100`
- source archive path: `proof_sources/mean_corrections/mean_corrections_body.tex`
- component SHA-256: `c494cebcaf476a01c5aae17d38f84dfc6b0092040c03b9aecbad1f5c3831ea29`
- associated checker: `proof_sources/mean_corrections/verify_exact.py`
- checker SHA-256 recorded in `navier_stokes_checks.json`: `5d54f54000e0cf03e70cfde1d298398bd29d091addf07887f70660af7c98a366`
- recorded checker result path: `proof_sources/mean_corrections/exact_checks.json`
- recorded result SHA-256: `5b81365d43a7108b16d11eba8cdf1dc267266226a0d90ff221ff67827f37495f`
- corrected 120-member bundle SHA-256: `43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5`
- DOI: `10.5281/zenodo.22678406`

The Kokuno repository did not expose a root `LICENSE` at the inspected public
revision. Therefore no checker implementation, proof body, or substantial
authored prose is copied. The target code is a clean-room implementation of the
short displayed operator identity only. Published source-side PASS metadata is
not treated as a local replay result.

## Replayed identity

For an exact local jet of `Psi_*`, the public MC14 formula gives

`Delta beta = -epsilon d_Z Psi_*`

and

`Delta gamma = (D_r + 1/R) Psi_*`.

The target replay verifies with exact `Fraction` arithmetic that the
cylindrical divergence contribution

`(D_r + 1/R) Delta beta + epsilon d_Z Delta gamma`

vanishes identically when mixed derivatives commute. The test fixture also
flips the required minus sign and obtains an exactly nonzero defect, so no
numerical tolerance can hide sign drift.

This is a local differential-operator regression aid for the existing Section
8 compact mean-repair lane. It does not construct the reserved bump, the
operator `T_1`, the full `A_1` remainder, the compact correction field, or an
actual oscillatory-wave debt.

## Actual target-side execution

Executed against the exact submitted source/test contents:

`PYTHONPATH=/mnt/data/kokuno7_mean_solenoidal/src python -m pytest -q -W error /mnt/data/kokuno7_mean_solenoidal/tests/test_kokuno_mean_solenoidal_replay.py`

Result: `5 passed in 0.06s`.

`python -W error -m py_compile /mnt/data/kokuno7_mean_solenoidal/src/openai_ns_reconstruction/kokuno_mean_solenoidal_replay.py /mnt/data/kokuno7_mean_solenoidal/tests/test_kokuno_mean_solenoidal_replay.py`

Result: exit `0`.

The exact fixture records:

- `Delta beta = -13/1700`
- `Delta gamma = -67/105`
- `D_r Delta beta = 19/2300`
- `d_Z Delta gamma = -371/1173`
- correct-sign divergence = `0`
- flipped-sign divergence = `-371/58650`
- tolerance = `0`

## Truth boundary

This finite replay does not establish compact-support existence, the imported
background/pulse hypotheses, the source checker's full result set, a materialized
Section 8 correction field, finite-stage residual improvement, all-stage
convergence, or the final force for one assembled witness.

`paper_exact_velocity_available=false`

`full_reconstruction=false`
