# Measured validation — 2026-09-10

## Scope and result

**91 tests passed locally**, with warnings treated as errors. The original 8 tests and the
4 provenance-manifest tests, 4 natural-axis tests and 3 background tests added on main were
preserved byte-for-byte. The natural-axis module and its provenance were merged from
main at `80d8282856d838c5fe0713411b2e54cb21fde2d3`; its composite-trapezoid formula was
retained with a NumPy 1.x/2.x-compatible implementation. The coefficient-radial-flux
formula and tests from `4f4abec6c5f63c53a62fc3949702fbcc62f5e883` were also integrated.
This validates the implemented components; **full reconstruction remains false**.

Environment: Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0. A wheel was built without network
access, installed into a separate environment, and imported from that wheel's site-packages
outside the source checkout. Available local NumPy/SciPy dependencies were shared explicitly;
this was not a fresh online dependency-resolution test. Wheel-based `demo` passed and
`audit --require-paper-exact` exited 2 as intended. The initial isolated smoke environment
could not find SciPy; after exposing the already-installed dependencies the wheel test passed.

Commands exercised:

```bash
python -m pytest -q -W error
python -m compileall -q src
python -m pip wheel . --no-deps --no-build-isolation -w dist
ns-reconstruct demo --output artifacts
ns-reconstruct audit --require-paper-exact   # expected exit code 2, not success
python examples/benchmark_quadrature.py
```

The CI workflow retains the four parallel diagnostic slices and adds full-suite coverage
for Python 3.10 / NumPy 1.x and Python 3.13 / NumPy 2.x, wheel installation outside the
checkout, and uploaded diagnostic artifacts. **Those matrix jobs are configured here;
local success is not a claim that remote CI has already passed.** Read the Actions run
for the actual commit's outcome.

## Numerical measurements

| Check | Measured value | Interpretation |
|---|---:|---|
| Manufactured coordinate q relative error, maximum in demo | 7.5162e-14 | Covers selected scales, not an all-input certificate |
| Toy singular-path log-log slope | -0.5049999999999998 | Matches leading exponent -0.505; profile remains toy |
| Exterior heat ODE defect, sampled maximum | 4.2588e-16 | Independently integrated H, H' and H'' |
| Exterior velocity/pressure NS residual norm | 6.0719e-8 | Finite-difference diagnostic on r>0, not the regular core |
| Independent Taylor–Green residual, step 0.04 | 3.57070e-4 | Analytically specified zero force |
| Same, step 0.02 | 8.92614e-5 | Approximately fourfold improvement |
| Same, step 0.01 | 2.23150e-5 | Consistent with second-order convergence on this test |

The demo writes exact runtime measurements, source pins, parameters, code hashes and
artifact hashes to `report.json`. Timing values vary by machine and run.

## Confirmed root-solver regression

At z=3e-9, t=0.9999999999999999, h=0.005, the old absolute-scale stopping condition returned
q=2.8532731732866517e-14. A 90-decimal-digit independent mpmath bisection gives
q=1.17258702425308461186340315510211794101754e-16. The revised solver returns
1.1725870242531042e-16, with relative error approximately 1.67e-14.
The old result was about 243 times the correct scale. The regression is not merely a
performance concern: q enters negative powers in the velocity.

Manufactured-root tests additionally cover q down to 1e-200. Direct-tau velocity tests avoid
the binary64 failure `1-tau == 1` for sufficiently small positive tau. Nonzero radii below
1e-15 are no longer silently forced to zero swirl.

## Radial-average microbenchmark

For U(X,eta)=eta exp(-X), X=2, eta=0.7, 1,000 calls per repeat and the median of five
repeats, the reproducible benchmark measured:

| Method | Median seconds | Absolute error |
|---|---:|---:|
| Previous 801-point composite trapezoid | 0.1186894 | 1.5762e-7 |
| Cached 32-point Gauss–Legendre rule | 0.0073503 | 5.5511e-17 |

The measured ratio was **16.15x for this single smooth radial-average operation only**.
It is not an end-to-end reconstruction speedup or a guarantee for arbitrary profiles.
Exact-average callbacks can avoid quadrature entirely. Changing the quadrature rule does
not establish regularity or support/cone/moment constraints for a profile.

## Correctness protections and unresolved work

Tests include pressure/viscosity signs, independent manufactured forces, bounded-time
stencils, invalid shapes/NaNs/step sizes, exact-axis conditions, potential-before-curl
localization, a deliberately non-axisymmetric swirl that has nonzero divergence,
cutoff short-circuit behavior, dyadic round trips and large exact integer covering powers.

Unresolved: complete regular leading profiles; recursive background coefficients;
transported phases and support separation; oscillatory stress and mean corrections;
convergent all-order iteration; paper-specific compact localization; smooth forcing
through t=1; full-field interval/formal certificates and local Lean builds.
No plot, fitted slope, residual-defined force or successful numerical test closes these gaps.
