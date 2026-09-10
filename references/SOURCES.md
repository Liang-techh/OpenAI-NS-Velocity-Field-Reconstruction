# Sources and provenance

Primary sources used for this reconstruction.

## OpenAI paper

**Finite Time Blowup for Navier–Stokes** (OpenAI, September 2026)

https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf

Key locations for the executable velocity reconstruction:

| Layer | Paper location |
|---|---|
| similarity variables | Eq. (4.1), Sec. 3.1/4.1 |
| leading velocity | Eq. (4.3) |
| Cartesian form | Eq. (4.5) |
| radial average | Eq. (4.6) |
| incompressibility / radial pressure | Eq. (4.7) |
| construction of leading profiles | Thm. 4.6 + App. A/B/C |
| all-order background coefficients | Eqs. (5.1)–(5.6) |
| vector-potential/streamfunction form | Eq. (5.27) |
| dyadic chart geometry | Eq. (6.1) and Sec. 6 |
| oscillatory stress realization | Sec. 7 |
| mean corrections | Sec. 8 |
| correction sequence | Eqs. (9.15)–(9.21) |
| summed local velocity | Eq. (9.21) |
| base vector potential | Eq. (10.1) |
| final localization | Eq. (10.4) |
| blow-up path | Eqs. (10.20)–(10.21) |

## Official Lean formalization

https://github.com/openai/NavierStokesAndEuler

The Lean repository is used as a theorem-level cross-check.  It is not itself a numerical implementation of `u(x,y,z,t)`.

## OpenAI announcement

https://openai.com/index/navier-stokes-solution/

The announcement gives the physical description: an inward-spiralling, axially stretched shrinking vortex whose velocity becomes unbounded while the external force remains smooth.

## Provenance policy

Every future concrete parameter/profile artifact should record:

1. paper equation/theorem/appendix location;
2. source commit of `openai/NavierStokesAndEuler` used for theorem cross-checking;
3. numerical choices required to instantiate existential constructions;
4. generated hashes for profile tables or coefficient files;
5. verification tolerances and truncation depth.

This prevents a numerically convenient surrogate from being mistaken for the published construction.
