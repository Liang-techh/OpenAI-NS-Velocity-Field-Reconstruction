# OpenAI NS Velocity Field Reconstruction

**Status: partial executable reconstruction, not the completed paper velocity field.**

The goal is an auditable numerical counterpart of the construction in OpenAI's
September 2026 paper *Finite Time Blowup for Navier–Stokes*. The repository
implements similarity geometry, leading-field kinematics, finite background
assembly, localization, pointwise moment-solver primitives, and numerical
verification. **The paper-specific leading profiles and the full correction
sequence are still missing.** The bundled Gaussian example is a **toy**, not
the OpenAI field. This is an independent project, not an OpenAI repository.

A passing test suite, a fitted exponent, or a force defined by a momentum
residual does **not** establish completion of the construction.

## Run from a clean checkout

Python 3.10 or newer is required. NumPy is the only runtime dependency.

```bash
python -m pip install -e ".[dev]"
python -m pytest -q -W error
python -m openai_ns_reconstruction status
python -m openai_ns_reconstruction verify --output verification.json
python -m openai_ns_reconstruction demo --output toy-demo --grid-size 41
```

`ns-reconstruct` is an equivalent installed command. Output paths must be new;
reports and demos do not silently overwrite existing files.

The demo writes `toy_velocity_slice.npz` and a JSON manifest with parameters,
a data SHA-256, source pins, and the explicit `toy-not-openai` label. On its
`y=0` slice, `velocity[j,i,k]` corresponds to `z[j]`, `x[i]`, and Cartesian
component `k`. Load it with `numpy.load(..., allow_pickle=False)`.

For an automated consumer that must refuse an incomplete reconstruction:

```bash
python -m openai_ns_reconstruction verify --require-paper-exact
```

**This intentionally returns exit code 2**, even when all numerical diagnostics
pass. Normal diagnostic failure returns 1; ordinary diagnostic success returns
0. Changing a profile's label does not bypass this completion gate.

## Implemented components

| Component | Implementation | Boundary |
|---|---|---|
| Similarity chart, Eqs. (4.1)-(4.2) | Log-scaled physical-branch root, direct `tau` interface, analytic chain rule | Float64 evaluation, not an interval certificate |
| Leading velocity and pressure, Eqs. (4.3)-(4.7) | Axis-regular evaluation, radial averages, pressure interface | Actual Theorem 4.6 profile constructor is missing |
| Radial quadrature | Cached Gauss–Legendre rule; optional analytic primitives | Resolution must be checked for each new profile |
| Background, Eqs. (5.1), (5.27) | Finite coefficient sum with analytic cutoff/curl product terms | Coefficients and admissible cutoff scales remain caller inputs |
| Appendix A primitives | Smooth-bump power moment matrix and zero-initialized quadratic moment iteration | Pointwise solve, not the full matching construction or uniform parameter bounds |
| Localization, Eq. (10.4) | `curl(c A) = c curl(A) + grad(c) cross A`, optional analytic callbacks | The direct swirl requires `c B` to be axisymmetric |
| Verification | Independent manufactured solutions, Taylor–Green refinement, divergence, pressure balance, finite-cylinder energy | No all-order smooth-force or global-energy proof |

Stages 3–6 (charts/phases, oscillatory realization, compact mean corrections,
and residual-improvement iteration) remain pending. See
[the stage ledger](docs/RECONSTRUCTION_PLAN.md),
[numerical details](docs/NUMERICS.md), and
[the next constructive tasks](docs/NEXT_TASKS.md).

## Near-singularity usage

Do **not** form `t = 1 - tau` when `tau` is below machine precision.
Use the explicit `*_from_tau` functions:

```python
from openai_ns_reconstruction import toy_gaussian_profile, blowup_probe

profile = toy_gaussian_profile()   # NOT OpenAI's constructed profile
v = blowup_probe(1e-80, 0.5, profile, h=0.005)
print(v.u_theta, v.q, v.X)
```

The geometry accepts `0 < h < 1/2`; the paper's narrower parameter domain is
`0 < h < 1/100`. The demo enforces the narrower domain. A representable
similarity coordinate does not guarantee that every derivative or velocity
fits in Float64; overflow and unrepresentable stencils are reported.

## Validation and reproducibility

`reports/diagnostics.json` records numerical errors, tolerances, refinement
steps and environment for this delivery. `reports/benchmark.json` compares
the original and new numerical methods on specified inputs. Timings are
local measurements, not general performance guarantees. Test logs and the
final delivery report record what was actually run.

To reproduce the benchmark, supply a local checkout of the pinned baseline:

```bash
python scripts/benchmark_baseline.py /path/to/baseline --output benchmark-new.json
```

The script checks the two baseline source blob hashes before comparison.
The CI configuration covers Python 3.10/3.11/3.13 and NumPy 1.x/2.x, but a
configuration is **not** evidence that remote CI has executed. This delivery
was locally tested on Python 3.13.5 with NumPy 2.3.5.

## Source of truth

- Paper: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- Official Lean repository: https://github.com/openai/NavierStokesAndEuler
- Announcement: https://openai.com/index/navier-stokes-solution/

Equation locations and provenance limits are recorded in
[references/SOURCES.md](references/SOURCES.md). The observed upstream commit
is metadata only: its Lean files have **not** been reviewed or compiled by
this delivery. No remote GitHub update is claimed.


### Concurrent changes retained

Before packaging, the delivery was integrated against remote snapshot
`b01aaeebc6ec8a39f6692109b6a3c81b92cd0f32` (five commits after the initial
benchmark baseline). Its source manifest, provenance tests and named
`slice`/`full` CI jobs are retained. The original and refreshed snapshots
were compared locally; the refreshed baseline's full Git tree hash was
verified against GitHub. Later remote changes still require normal review.
