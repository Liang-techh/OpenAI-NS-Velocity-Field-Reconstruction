# Endpoint residual-jet limit majorant provenance

Status: **formal-structure only**. This document does not claim a paper-exact residual, endpoint jet, forcing, or velocity field.

## Upstream source boundary

The repository pins OpenAI `NavierStokesAndEuler` at commit `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`. In `NavierStokes/CandidateFromLimits.lean`, the final Section 10 construction assumes locally uniform `t -> 1-` limits for every full spacetime derivative of the closed-past Navier--Stokes residual. Existing executable modules already consume those limit tensors (`endpoint_jets.py`, `traced_residual.py`) and build the conditional Spatial-Borel scale/tail machinery (`endpoint_scale_schedule.py`), but none of those modules proves that the actual residual has the required limits.

`src/openai_ns_reconstruction/endpoint_limit_majorant.py` adds a reusable sufficient-condition layer for that missing hypothesis. It is not a replacement theorem from the Lean source; it is an elementary fundamental-theorem-of-calculus implication that can later discharge a `CandidateFromLimits` limit obligation once the real localized residual is available.

For one derivative degree `n` and one compact spatial window `K_m`, the caller must independently prove

`sup_{x in K_m} ||partial_t D^n R(t,x)|| <= C (T-t)^(-alpha)`

for `valid_from <= t < T`, with `C >= 0` and `0 <= alpha < 1`. The module then computes the uniform Cauchy bound

`C/(1-alpha) * ((T-t)^(1-alpha) - (T-s)^(1-alpha))`

for `t <= s < T`, together with the endpoint tail modulus

`C/(1-alpha) * (T-t)^(1-alpha)`.

The default `valid_from=3/4` matches the already-landed Section 10 time switch regime where the executable switch is identically one, but callers may provide another rigorously justified interval.

## Independent checks

`tests/test_endpoint_limit_majorant.py` uses a labelled manufactured family with an analytically known endpoint limit. Its spatial amplitude attains the supplied compact-window supremum, so the endpoint bound is sharp for that fixture. A separate midpoint quadrature of `C(1-t)^(-alpha)` cross-checks the closed-form interval integral. Invalid nonintegrable exponents (`alpha >= 1`), negative coefficients, out-of-regime times, and malformed derivative/window indices fail closed.

The manufactured family is a regression oracle only. It is not an OpenAI velocity/residual surrogate and is never passed to the paper-exact gate.

## What remains unresolved

This increment does **not** derive `C` or `alpha` from the actual localized field, does not infer them from finite samples or fitted slopes, and does not produce any endpoint tensor by itself. The full construction still needs, for every derivative degree and every compact window, genuine analytic estimates for the actual closed-past residual that satisfy a limit criterion of this kind (or an equivalent proof of the `CandidateFromLimits` hypotheses). Genuine all-order template derivative bounds are also still required for the Spatial-Borel scale schedule and smooth force extension through `t=1`.

Accordingly Stage 7 remains `formal-structure` and `paper_exact_velocity_available` must remain `false`.
