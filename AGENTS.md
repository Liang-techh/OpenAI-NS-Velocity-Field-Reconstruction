# Contributor / agent contract

The goal is the full paper construction, not a visually plausible blow-up surrogate.
Read README.md, docs/RECONSTRUCTION_PLAN.md and references/provenance_manifest.json first.

- Preserve toy/formula/structure/diagnostic distinctions. Never promote status from a green test alone.
- Every new mathematical constructor needs paper locations, parameter choices, domain/support conditions and independent tests.
- Do not call `f=R(u,p); R-f=0` independent verification, or infer smooth forcing from samples.
- Do not silently set missing pressure/coefficients to zero, regularize the exterior at the axis, or claim all finite coefficients imply convergence.
- Keep the direct swirl and localization axisymmetric when invoking divergence-free composition.
- Use direct tau APIs near t=1. Convergence tests must vary quadrature, spatial/time steps and truncation independently.
- Run `python -m pytest -q -W error`, `ns-reconstruct demo --output artifacts`, and confirm `ns-reconstruct audit --require-paper-exact` exits 2 while constructors are incomplete.
- Preserve existing tests and reconcile source/manifest/CI changes against the latest main; never force-push over another agent.
- Record measured results and remaining blockers. Do not invent PDF hashes, Lean builds or CI outcomes.
