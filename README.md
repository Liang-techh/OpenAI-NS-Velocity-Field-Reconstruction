# OpenAI NS Velocity Field Reconstruction

An independent, **partial** executable reconstruction of the velocity construction in
OpenAI's September 2026 *Finite Time Blowup for Navier–Stokes*. This is not an OpenAI repository.

**Current result:** tested similarity kinematics; natural-axis/range structure, an admissible
ideal-prefix pressure witness, the constructive outgoing scalar schedule, executable
`finalAngular/clockWeight`, the actual outgoing `SchedulePressure.axisPressure`, an analytic
no-sampling low-|Z| `H^2` margin instantiating the theorem-side `sigma`, the theorem-faithful
`Lambda/C` selection algebra, a fail-closed AxisContraction-style norm-ledger adapter, and the
canonical analytic-input propagation that turns an independently certified common complex radius
`rho` and common value bound `B` into `epsilon=rho/2`, `radiusLoss(1/2)=12`, coefficient norm
`<=12B`, normalized amplitude `M=12`, and a conservative complete naturalResolvent factorial-series
bound; Section 5 pressure/Omega recurrence rows, the Eq. (5.7) singular Picard primitive, Lemma 5.2
compact five-moment repair, Eq. (5.15) forward reconstruction from supplied repaired analytic
coefficient data, a finite-prefix SlowBorel/DiagonalScale cutoff-scale constructor, exact finite-
prefix plateau/transition/zero-tail support certification, and a theorem-shaped fixed-prefix
truncation/tail-order arithmetic gate; fixed dyadic geometry, the slow-label/base-jet bridge,
executable squared slow partitions, a constructive auxiliary-slot centers/`r0` witness, the
physical-support-to-`SlotColoring.Adj` bridge, pointwise and uniform phase-normal adapters, and the
pinned BasePhaseGeometry frame/damping implication bounds with the corrected family-level
`phaseConstant(M)=normalConstant(frequencyBound(M))`; spatial/time localization, the pinned closed-
past `zeroBefore` branch, endpoint Taylor–Borel right-extension infrastructure, the full-spacetime-
to-normal endpoint-jet adapter, SpatialBorel scale arithmetic, a value-level CandidateFromLimits
traced-residual/Borel glue bridge, and the exact Section 10 support-cylinder/fixed-time energy
implication; plus reproducible numerical diagnostics.

**Not yet delivered:** an actual admissible common complex analytic neighborhood `rho`, the true
common complex-field sup bound `B`, and complex `realPartSup` for the paper's coefficient family;
the coefficient-space fixed-point fields and complete regular leading profile; the profile-derived
converged all-order background solve, true eta-dependent repaired hierarchy/support closure,
uniform `C[j,m]` bounds, finite residual-cancellation identities, and completed infinite recursive
cutoff/all-jets-flat argument; the paper-exact Proposition 5.5 base-field provider and actual
LocalBase/vector hypotheses on every active box, stress-cone/amplitude/wave construction, mean
corrections and convergent correction sequence; or the final compact field with actual closed-past
residual full-spacetime derivative limits, genuine analytic template majorants, a force proved
smooth through the singular time, and the paper's uniform bounded-energy conclusion along the
actual blow-up limit.

A successful demo or a green test suite is **not** a reconstruction of the full counterexample.
The CLI reports `full_reconstruction: false` / `paper_exact_velocity_available: false` and refuses
`--require-paper-exact` requests.

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
The `status` and `audit` commands derive their stage truth from the same fail-closed runtime source,
while retaining their respective output schemas. The demo produces labelled diagnostic/toy
artifacts only; those outputs are not samples of the final counterexample.

To demand a complete paper reconstruction:

```bash
ns-reconstruct audit --require-paper-exact
```

This currently exits **2**, intentionally. It is a conservative completion gate, not a formal
proof checker. `demo --require-paper-exact` also refuses before generating toy data.

## Implemented components

| Component | Implementation and boundary |
|---|---|
| Similarity coordinates, Eq. (4.1) | Relative-scale root solve, finite-input checks, direct `tau` API |
| Leading velocity, Eqs. (4.3)–(4.7) | Caller-supplied profiles, regular-axis handling, pressure evaluation; natural-axis/range formulas, outgoing schedule, actual `axisPressure`, analytic low-|Z| margin, theorem-faithful `Lambda/C` selection, theorem-shaped conservative `remainderBound/remainderLip` propagation, and canonical analytic-input norm/resolvent propagation are present. Given independently certified `rho,B`, the landed adapter fixes `M=12`, common fixed-field norm `<=12B`, and the complete factorial-series majorant. **The actual admissible `rho`, true common complex-field bound `B`, `realPartSup`, coefficient-space fixed point, completed `NaturalProfileAssembly`, and final support/moment/matching/cone certificates are still missing.** |
| Profile averages | Cached 32-point Gauss–Legendre rule; optional exact-average callbacks |
| Exterior heat swirl, Appendix A.6 | Adaptive evaluation of H and derivatives, swirl and centrifugal pressure; **r>0 only** |
| Section 5 background rows | Eq. (5.2) radial flux, Eq. (5.27) streamfunction/vector potential, Eq. (5.5) pressure row, regular Eq. (5.6) `Omega_k/X`, Eq. (5.7) singular inverse/Picard primitive, Lemma 5.2 compact five-moment repair, Eq. (5.15) forward `F_n/V_n/Pi_n` reconstruction from supplied repaired analytic data, finite-prefix recursive cutoff scheduling from supplied analytic `C[j,m]`, exact finite-prefix cutoff support/truncation certification, and exact fixed-prefix tail-order arithmetic for `2^-J q^(h(J+1)-m)` / physical exponent `h(J+1)+b-2M`. **The paper-derived `A0/A1/f_n`, converged Lemma 5.1 coefficient hierarchy, true eta-dependent support/stress closure, true uniform `C[j,m]`, completed infinite schedule/theorem-level local-finiteness argument, order-by-order residual cancellation identities, and Proposition 5.3 all-jets-flat decay are not yet constructed.** |
| Dyadic geometry / phase, Sections 6–7 | Fixed chart scaling and exact covering matrices; active-shell/slow-label bridge; squared partitions; constructive 2250-color/rational-center/common-`r0` witness; physical slow-support-to-`SlotColoring.Adj` bridge with cross-band common-point handling; Section 7.1 phase/tangent-frame algebra; pointwise rounded-normal adapters; a fail-closed UniformLocalBase implication bridge; and pinned BasePhaseGeometry frame/damping consequence bounds, including the enlarged family-level phase constant. **The paper-exact Proposition 5.5 provider and actual LocalBase C1/C2, unit/orthogonality, normal-closeness and slot-normal derivative-closeness certificates remain missing; the landed implications cannot be promoted until those true hypotheses are supplied, after which stress/amplitude/wave work still remains.** |
| Section 10 localization | Official spatial/time support/plateau geometry represented by explicit C-infinity bumps, analytic spatial cutoff gradient, time-switch derivative, activated-field adapters and independent residual-identity cross-check; pinned closed-past `zeroBefore/pastVelocity/pastPressure` branch with diagnostic residual evaluation; endpoint Taylor–Borel infrastructure, full-spacetime endpoint-jet adapter, exact-rational SpatialBorel scale scheduling from supplied analytic derivative bounds, and a value-level traced-residual/Borel glue bridge for supplied endpoint tensors. The support certificate also records the exact radius-`1/4`, height-`1/2` cylinder of volume `pi/32`, hence `E(t)<=pi M^2/64` from an independently certified fixed-time `|u|<=M`. **Transition-collar point values are not claimed equal to Mathlib's noncomputable bump; actual closed-past residual derivative limits, true analytic template majorants, smooth force gluing through `t=1`, and a uniform bounded-energy estimate along the actual blow-up limit are still missing.** |
| Verification | Independent manufactured-solution and refinement checks; bounded-time finite differences |
| Audit/provenance | Source-pinned ledger, explicit blockers, synchronized `status`/`audit` truth surface, fail-closed completion gate |

### Near-singularity evaluation

Use `tau=1-t` directly when working close to the singular time. Forming `t=1-tau` in binary64
can round to 1 and irreversibly lose small positive tau.

```python
from openai_ns_reconstruction import similarity_coordinates_from_tau
from openai_ns_reconstruction.heat_exterior import HeatExterior

s = similarity_coordinates_from_tau(r=1e-50, z=0.0, tau=1e-100, h=0.005)
print(s.q, s.X)

exterior = HeatExterior(h=0.005, c_inf=1.0)
print(exterior.swirl_from_tau(r=1.0, tau=0.2))
print(exterior.pressure_from_tau(r=1.0, tau=0.2))
```

`HeatExterior` implements a parameterized published exterior component, **not** an admissible
smooth-axis `LeadingProfile` and not the full field. The numerical code accepts `0<h<1/2` where
the coordinate formulas make sense; the paper's complete construction imposes the stricter
`0<h<1/100` and additional parameter conditions.

## Verification boundaries

Finite differences, quadrature, and sampled convergence are diagnostics, not interval bounds or
Lean certificates. Setting `f=R(u,p)` and checking `R-f` with the same stencil is tautological;
our regression suite additionally uses independently specified manufactured forces. A generic
`B(x,y,z,t)e_theta` need not be divergence-free: the direct swirl must be axisymmetric, and
localization must preserve that property. Prefer the `LocalField.from_axisymmetric` and
`LocalizedField.from_axisymmetric` factories.

The regular inner core and materialized coefficient-space fixed point, actual analytic-neighborhood
and complex-field sup witnesses, `realPartSup`, profile-derived Section 5 source and converged
hierarchy, true uniform coefficient-template bounds and completed infinite recursive cutoff
schedule/residual cancellation, paper-exact base-field/slow-box and vector hypotheses,
stress/wave realization, mean correction, infinite-order summation, completed compact field,
actual closed-past residual full spacetime jets and locally uniform endpoint derivative limits,
genuine analytic compact-template derivative majorants, smooth-force extension, blow-up-path
closure, and uniform bounded-energy conclusion remain explicit blockers. The fixed-time
compact-support energy implication does not replace that last uniform estimate.
See [the reconstruction plan](docs/RECONSTRUCTION_PLAN.md),
[the measured validation report](docs/VALIDATION_2026-09-10.md), and
[the machine-readable manifest](references/provenance_manifest.json).

## Sources

- Paper: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- Official Lean source: https://github.com/openai/NavierStokesAndEuler
- Pinned upstream commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Announcement: https://openai.com/index/navier-stokes-solution/

Pinning the source commit does not mean its Lean code has been built or its theorem-to-runtime
mapping verified here. The paper was inspected in rendered pages; its binary hash has not been
computed. These limitations are recorded rather than filled with guessed provenance.
