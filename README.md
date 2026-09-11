# OpenAI NS Velocity Field Reconstruction

An independent, **partial** executable reconstruction of the velocity construction in
OpenAI's September 2026 *Finite Time Blowup for Navier–Stokes*. This is not an OpenAI repository.

**Current result:** tested similarity kinematics; natural-axis/range structure and an admissible
ideal-prefix pressure witness; Section 5 pressure/Omega recurrence rows; fixed dyadic chart and
pointwise phase-estimate adapters; spatial/time localization structure; and reproducible numerical
diagnostics.
**Not yet delivered:** the complete regular leading profile, the coupled all-order background solver,
stress-cone/amplitude/wave construction, mean corrections, convergent correction sequence, or the
final compact field with a force proved smooth through the singular time.

A successful demo or a green test suite is **not** a reconstruction of the full counterexample.
The CLI reports `full_reconstruction: false` and refuses `--require-paper-exact` requests.

## Run from a clean checkout

Requires Python 3.10+; NumPy and SciPy are installed as dependencies.

```bash
python -m pip install -e '.[dev]'
python -m pytest -q -W error
ns-reconstruct audit
ns-reconstruct demo --output artifacts
```

`python -m openai_ns_reconstruction` is equivalent to the console command.
The demo produces `report.json`, `toy_blowup_probe.csv`, `toy_velocity_samples.csv`,
and `heat_exterior.csv`. The report contains parameters, package/runtime versions,
source pins, numerical tolerances, and SHA-256 hashes for code and generated data.
It runs offline after installation. CSV names explicitly distinguish toy samples from
parameterized paper components. They are not samples of the final counterexample.

To demand a complete paper reconstruction:

```bash
ns-reconstruct audit --require-paper-exact
```

This currently exits **2**, intentionally. It is a conservative completion gate, not a
formal proof checker. `demo --require-paper-exact` also refuses before generating toy data.

## Implemented components

| Component | Implementation and boundary |
|---|---|
| Similarity coordinates, Eq. (4.1) | Relative-scale root solve, finite-input checks, direct `tau` API |
| Leading velocity, Eqs. (4.3)–(4.7) | Caller-supplied profiles, regular-axis handling, pressure evaluation; natural-axis/range formulas and an abstract admissible pressure witness are executable, but the paper's full profile/schedule/fixed-point data are still missing |
| Profile averages | Cached 32-point Gauss–Legendre rule; optional exact-average callbacks |
| Exterior heat swirl, Appendix A.6 | Adaptive evaluation of H and derivatives, swirl and centrifugal pressure; **r>0 only** |
| Section 5 background rows | Eq. (5.2) radial flux, Eq. (5.27) streamfunction/vector potential, Eq. (5.5) pressure row, and regular Eq. (5.6) `Omega_k/X` from supplied coefficient jets; **no coupled (5.3)–(5.4) solver, Lemma 5.2 repair, or paper cutoff schedule yet** |
| Dyadic geometry / phase, Sections 6–7 | Fixed chart scaling and exact covering matrices; Section 7.1 phase/tangent-frame algebra, scalar PhaseEstimates gate, and pointwise `rounded_normal_estimates` adapter; **no uniform slow-box certificate or transported-wave construction** |
| Section 10 localization | Official spatial/time support/plateau geometry represented by explicit C-infinity bumps, analytic spatial cutoff gradient, time-switch derivative, activated-field adapters, and independent numerical residual-identity cross-check; **transition-collar point values are not claimed equal to Mathlib's noncomputable bump and the smooth force extension through `t=1` is missing** |
| Verification | Independent manufactured-solution and refinement checks; bounded-time finite differences |
| Audit/provenance | Source-pinned ledger, explicit blockers, fail-closed completion gate |

### Near-singularity evaluation

Use `tau=1-t` directly when working close to the singular time. Forming `t=1-tau` in
binary64 can round to 1 and irreversibly lose small positive tau.

```python
from openai_ns_reconstruction import similarity_coordinates_from_tau
from openai_ns_reconstruction.heat_exterior import HeatExterior

s = similarity_coordinates_from_tau(r=1e-50, z=0.0, tau=1e-100, h=0.005)
print(s.q, s.X)

exterior = HeatExterior(h=0.005, c_inf=1.0)
print(exterior.swirl_from_tau(r=1.0, tau=0.2))
print(exterior.pressure_from_tau(r=1.0, tau=0.2))
```

`HeatExterior` implements a parameterized published exterior component, **not** an
admissible smooth-axis `LeadingProfile` and not the full field. The numerical code accepts
`0<h<1/2` where the coordinate formulas make sense; the paper's complete construction
imposes the stricter `0<h<1/100` and additional parameter conditions.

## Verification boundaries

Finite differences, quadrature, and sampled convergence are diagnostics, not interval bounds
or Lean certificates. Setting `f=R(u,p)` and checking `R-f` with the same stencil is tautological;
our regression suite additionally uses independently specified manufactured forces.
A generic `B(x,y,z,t)e_theta` need not be divergence-free: the direct swirl must be
axisymmetric, and localization must preserve that property. Prefer the
`LocalField.from_axisymmetric` and `LocalizedField.from_axisymmetric` factories.

The regular inner core, profile matching/moments/cone constraints, coupled coefficient solver,
compact moment repair and recursive cutoff schedule, uniform slow-box estimates, stress/wave
realization, mean correction, infinite-order summation, completed compact field, and smooth-force
extension remain explicit blockers. See [the reconstruction plan](docs/RECONSTRUCTION_PLAN.md),
[the measured validation report](docs/VALIDATION_2026-09-10.md), and
[the machine-readable manifest](references/provenance_manifest.json).

## Sources

- Paper: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- Official Lean source: https://github.com/openai/NavierStokesAndEuler
- Pinned upstream commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Announcement: https://openai.com/index/navier-stokes-solution/

Pinning the source commit does not mean its Lean code has been built or its theorem-to-runtime
mapping verified here. The paper was inspected in rendered pages; its binary hash has not
been computed. These limitations are recorded rather than filled with guessed provenance.
