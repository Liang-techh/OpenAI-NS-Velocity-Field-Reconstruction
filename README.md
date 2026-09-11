# OpenAI NS Velocity Field Reconstruction

An independent, **partial** executable reconstruction of the velocity construction in
OpenAI's September 2026 *Finite Time Blowup for Navier–Stokes*. This is not an OpenAI repository.

**Current result:** tested similarity kinematics; natural-axis/range structure, an admissible
ideal-prefix pressure witness, the constructive outgoing scalar schedule, executable
`finalAngular/clockWeight`, the actual outgoing `SchedulePressure.axisPressure`, and an analytic
no-sampling low-|Z| `H^2` margin instantiating the theorem-side `sigma`; Section 5
pressure/Omega recurrence rows, the Eq. (5.7) singular Picard primitive, the Lemma 5.2 compact
five-moment repair, and a finite-prefix SlowBorel/DiagonalScale cutoff-scale constructor from
supplied analytic coefficient bounds; fixed dyadic geometry, the slow-label/base-jet bridge,
executable squared slow partitions, a constructive auxiliary-slot centers/`r0` witness under the
discrete interaction hypotheses, and pointwise phase-estimate adapters; spatial/time
localization, endpoint Taylor–Borel right-extension infrastructure, the full-spacetime-to-normal
endpoint-jet adapter, and the SpatialBorel boundSum/localScale/doubling schedule arithmetic from
supplied analytic template bounds; plus reproducible numerical diagnostics.
**Not yet delivered:** the complete regular leading profile, constructive `Lambda/C` choices and
coefficient-space fixed-point fields; the profile-derived converged all-order background solve,
full Eq. (5.15) repaired-coefficient reconstruction, true uniform `C[j,m]` bounds and completed
infinite recursive cutoff/residual-decay argument; the physical-support-to-`SlotColoring.Adj`
bridge needed to turn the landed slot witness into full Lemma 6.1 separation, the paper-exact
base-field provider and uniform phase estimates, stress-cone/amplitude/wave construction, mean
corrections and convergent correction sequence; or the final compact field with actual residual
endpoint derivative limits, genuine analytic template majorants, and a force proved smooth
through the singular time.

A successful demo or a green test suite is **not** a reconstruction of the full counterexample.
The CLI reports `full_reconstruction: false` / `paper_exact_velocity_available: false` and
refuses `--require-paper-exact` requests.

## Run from a clean checkout

Requires Python 3.10+; NumPy and SciPy are installed as dependencies.

```bash
python -m pip install -e '.[dev]'
python -m pytest -q -W error
ns-reconstruct status
ns-reconstruct audit
ns-reconstruct demo --output artifacts
```

`python -m openai_ns_reconstruction` is equivalent to the console command.
The `status` and `audit` commands derive their stage truth from the same fail-closed runtime
source, while retaining their respective output schemas. The demo produces
`report.json`, `toy_blowup_probe.csv`, `toy_velocity_samples.csv`, and `heat_exterior.csv`.
The report contains parameters, package/runtime versions, source pins, numerical tolerances,
and SHA-256 hashes for code and generated data. It runs offline after installation. CSV names
explicitly distinguish toy samples from parameterized paper components. They are not samples
of the final counterexample.

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
| Leading velocity, Eqs. (4.3)–(4.7) | Caller-supplied profiles, regular-axis handling, pressure evaluation; natural-axis/range formulas, an abstract admissible pressure witness, the pinned Gaussian-flat outgoing schedule (`S=32` derivative-bound witness, `flattenLength`, `shapeExponent`), complete executable `finalAngular/clockWeight`, the actual all-real-line `SchedulePressure.axisPressure`, and an analytic low-|Z| `H^2` margin yielding the theorem-side `sigma=sqrt(m)/20` are present; **`Lambda/C`, the coefficient-space fixed-point fields, completed `NaturalProfileAssembly`, and final support/moment/matching/cone certificates are still missing** |
| Profile averages | Cached 32-point Gauss–Legendre rule; optional exact-average callbacks |
| Exterior heat swirl, Appendix A.6 | Adaptive evaluation of H and derivatives, swirl and centrifugal pressure; **r>0 only** |
| Section 5 background rows | Eq. (5.2) radial flux, Eq. (5.27) streamfunction/vector potential, Eq. (5.5) pressure row, regular Eq. (5.6) `Omega_k/X`, the Eq. (5.7) six-component singular inverse/Picard-map primitive, the Lemma 5.2 / Eqs. (5.14)–(5.16) compact five-moment repair with two U bumps and three E bumps, and a finite-prefix `SlowBorelBase`/`DiagonalScale` cutoff schedule from supplied analytic normalized-template bounds `C[j,m]`; **the paper-derived `A0/A1/f_n`, converged Lemma 5.1 coefficient, full Eq. (5.15) repaired `F_n/V_n/Pi_n` reconstruction, true uniform `C[j,m]`, completed infinite schedule/local-finiteness argument, and arbitrary-order residual decay are not yet constructed** |
| Dyadic geometry / phase, Sections 6–7 | Fixed chart scaling and exact covering matrices; Eq. (6.8) active-shell/slow-label bridge with typed `TangentialBaseJetProvider`; paper-admissible normalized C-infinity translate squared partitions for dyadic q and the `S_*^-3` `(R,Z,T)` slow mesh; a constructive Section 6.2 2250-color/rational-center/common-`r0` slot witness conditional on the discrete interaction hypotheses; Section 7.1 phase/tangent-frame algebra, scalar PhaseEstimates gate, and pointwise `rounded_normal_estimates` adapter; **the physical-support/enlargement → `SlotColoring.Adj` bridge and hence full Lemma 6.1 separation, paper-exact provider instantiation, and uniform slow-box certificate are not yet complete** |
| Section 10 localization | Official spatial/time support/plateau geometry represented by explicit C-infinity bumps, analytic spatial cutoff gradient, time-switch derivative, activated-field adapters and independent numerical residual-identity cross-check; endpoint Taylor–Borel infrastructure reproduces `doublingEnvelope`, locally finite future evaluation, endpoint value and `t>=T+1` zero support; a dense full-spacetime endpoint-jet adapter contracts every derivative slot with the pinned time direction before Borel evaluation; an exact-rational SpatialBorel `boundSum/localScale` schedule from supplied analytic compact-template derivative bounds certifies the corresponding `2^-j` tail inequalities; **transition-collar point values are not claimed equal to Mathlib's noncomputable bump, actual residual endpoint derivative limits and true analytic template majorants are missing, and smooth force gluing through `t=1` is not proved** |
| Verification | Independent manufactured-solution and refinement checks; bounded-time finite differences |
| Audit/provenance | Source-pinned ledger, explicit blockers, synchronized `status`/`audit` truth surface, fail-closed completion gate |

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

The regular inner core, `Lambda/C` and fixed-point profile data, profile-derived Section 5
matrices/source and converged coefficient solve, full Eq. (5.15) reconstruction after the
landed compact five-moment repair, true uniform coefficient-template bounds and completed
infinite recursive cutoff schedule, the physical-support-to-adjacency bridge completing Lemma
6.1 after the landed squared partitions/slot witness, uniform slow-box estimates, stress/wave
realization, mean correction, infinite-order summation, completed compact field, actual
residual endpoint full spacetime jets and locally uniform derivative limits, genuine analytic
compact-template derivative majorants, and smooth-force extension remain explicit blockers.
See [the reconstruction plan](docs/RECONSTRUCTION_PLAN.md),
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
