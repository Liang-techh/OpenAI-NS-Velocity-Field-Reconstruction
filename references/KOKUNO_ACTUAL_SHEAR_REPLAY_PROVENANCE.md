# Kokuno corrected actual-shear replay provenance

## Scope

This record accompanies one independently authored exact checker for the short
corrected-edition identity

`-g·T = -|g0| T_N - (g-g0)·T`

in the oscillatory pulse calculation.  It verifies only the Euclidean
coordinate split in the public `(N,K)` shear plane and makes omission of the
actual-shear remainder fail closed when that remainder is nonzero.

It does **not** copy or execute the Kokuno checker, construct the oscillatory
velocity, establish the covariance/curl/shear estimates, prove positive
production, solve the mean-correction equations, validate the finite-stage
cycle, or certify any all-stage convergence claim.

## Public source identity

- Source repository: `KokunoYumeto/yang-mills-interacting-workbench`
- Public workbench commit inspected: `fab69fdc4ac197159b8e6ae8d73a82bde2b20d55`
- `navier-stokes/navier_stokes_workbench.tex` Git blob SHA-1:
  `205a99807302e21a51c5eaf223390c0dfc42bcd0`
- Frozen research-state source commit:
  `e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f`
- Corrected release DOI: `10.5281/zenodo.22678406`
- Corrected 120-member bundle SHA-256 (published metadata):
  `43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5`
- `NS-oscillations` body SHA-256:
  `3fd5c61de39b6f65a581ead5b29c741c2e0f46fb30f62e582dc24d93bbe83615`
- Recorded source checker `proof_sources/oscillations/exact_checks.py` SHA-256:
  `a44ba49990c3b604c883d5b4e1726cb33068867cd85b23f3a2c22847e4ce8184`
- Recorded checker result `proof_sources/oscillations/exact_checks_results.json` SHA-256:
  `30d7a4d934e984a77cfcb04188d72c35391c46e38bcbf25a2024e38bcb783a12`
- Source manuscript capture SHA-256:
  `8c8a94ad9ac824c8b605b9827cadf7beaca48bd10b380de3cfc872a2c37afa81`
- Source pages recorded for `NS-oscillations`: `62–87,157–165`.

The bundle/checker/result hashes above are provenance metadata.  The binary ZIP
could not be decoded by the GitHub text connector in this run, so they are not
reported as fresh local byte re-hashes and the original checker was not
executed.

## Copyright / license boundary

No clear repository-level reuse license was found for the Kokuno workbench.
Accordingly this increment copies no checker implementation and no substantial
source prose.  The Python code is an independent implementation of the short
public mathematical identity and exact rational fail-closed behavior.

## Independent replay semantics

Write

- `g0 = m N`, with exact `m=|g0| >= 0`;
- `g-g0 = d_N N + d_K K`;
- `T = T_N N + T_K K`.

Using exact `Fraction` arithmetic, the checker compares

- actual production: `-((m+d_N) T_N + d_K T_K)`, and
- corrected split: `-m T_N -(d_N T_N + d_K T_K)`.

Their exact residual must be zero.  The difference between actual production
and the frozen leading term is exactly the retained remainder.  The checker
therefore rejects a leading-only substitution whenever that exact remainder is
nonzero; it has no numerical tolerance or near-zero gate.

## Truth boundary

This replay is one algebraic regression aid for the corrected actual-shear
term.  It is not the source `NS-oscillations` checker, not independent
validation of all oscillatory identities, not completion of target task NS029,
and not evidence for NS030/NS031 or all-stage Section 9 convergence.

`paper_exact_velocity_available=false`

`full_reconstruction=false`
