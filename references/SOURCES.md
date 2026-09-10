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


## Source snapshot used by the 0.2.0 delivery

- User repository baseline: `b0d964ab4a69c25863078830faf76296550ae4ea`.
- Upstream commit metadata observed: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
- Paper pages inspected for this change include printed pages 24–25
  (Eqs. 4.1–4.7), 45 (Eq. 5.1), 56 (Eq. 5.27) and 127 (Lemma A.2).
- No paper-byte checksum was established. No Lean source review or Lean build
  was completed; the observed upstream SHA is **not** a theorem-level validation.
- Numerical measurements in `reports/` are local reproduction evidence only.
- New C-infinity cutoffs and the Gaussian test profile are independent numerical
  choices, not the paper's complete instantiated parameter/profile data.

- Integration base refreshed to `b01aaeebc6ec8a39f6692109b6a3c81b92cd0f32`;
  reconstructed Git tree `8788a0faf709c02ece63a41774d29f9754d7a59b` was checked
  byte-for-byte against the remote tree. The benchmark still intentionally
  compares against the original `b0d964a` source files.
