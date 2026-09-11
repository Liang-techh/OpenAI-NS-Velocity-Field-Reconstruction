"""Machine-readable scope and source pins; passing tests cannot change truth status.

This runtime status is intentionally conservative. It mirrors the repository's
current implemented *boundaries* closely enough for the CLI to stay useful, but it
never upgrades a construction stage to paper-exact merely because code or tests
exist. The detailed source-of-truth ledger remains
``references/provenance_manifest.json``.
"""
from __future__ import annotations
import copy

SOURCE_PINS = {
    "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction",
    # Stable historical integration facts, not a claim that these are current main.
    "pr5_integration_merge_commit": "c0a08e68b0f65a59a4dca70c50934957a0608475",
    "integrated_pr5_head": "9c7ad1ff6dd4b21806ff5ee2248974975542b76d",
    "paper_url": "https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf",
    "paper_equations_checked": [
        "4.1-4.7",
        "5.1-5.6",
        "5.27",
        "6.1-6.6",
        "7.1-7.11",
        "10.4",
        "Lemma 5.1",
        "Lemma A.2",
    ],
    "upstream_lean_repository": "openai/NavierStokesAndEuler",
    "upstream_commit_metadata_observed": "f9e8bc5b38b6e212696e8a30e3e91517af887bbd",
    "lean_source_reviewed": True,
    "lean_compiled": False,
    "paper_hash_verified": False,
}

STAGES = [
    {"id": 0, "name": "similarity geometry", "status": "paper-exact-formula",
     "implemented": "Eq. (4.1) similarity chart with stable direct-tau evaluation and analytic coordinate derivatives.",
     "remaining": "No interval or arbitrary-precision runtime certificate."},
    {"id": 1, "name": "leading profile", "status": "formal-structure",
     "implemented": (
         "Leading kinematics, natural-axis polynomials/rescaling, heat exterior, pointwise moment primitives, "
         "NaturalAxisRange h/j bounds, unique H-root bracket, delta=j/10 and conditional sigma=sqrt(m)/20."
     ),
     "remaining": (
         "Materialize the coefficient-space fixed-point phi/u/average/pressure fields; construct the admissible "
         "pressure schedule and uniform H^2 margin; instantiate Lambda/C; then verify support, moments, matching and cone conditions."
     )},
    {"id": 2, "name": "all-order background", "status": "formal-structure",
     "implemented": (
         "Finite coefficient assembly, Eq. (5.2) radial flux, Eq. (5.27) streamfunction/vector potential, "
         "analytic cutoff/curl terms, and the Eq. (5.5) pressure recurrence row with Pi_n(0,eta)=0 integration."
     ),
     "remaining": (
         "Coupled Eqs. (5.3)-(5.4) recursive solve, Eq. (5.6) Omega source construction, Lemma 5.2 compact moment repair, "
         "paper-selected recursive cutoff schedule and all-order residual estimates."
     )},
    {"id": 3, "name": "dyadic charts and transported phases", "status": "formal-structure",
     "implemented": (
         "Fixed dyadic geometry, Section 7.1 phase/wavevector/tangent-frame algebra, and the pinned PhaseEstimates "
         "scalar scale gate for the Eq. (7.9)-(7.11) error envelope."
     ),
     "remaining": (
         "Instantiate paper-exact base fields/slow supports and certify local C2/base-field hypotheses so the scalar gate "
         "becomes a uniform active-box phase/frame/damping certificate."
     )},
    {"id": 4, "name": "oscillatory stress realization", "status": "pending",
     "remaining": "Stress-cone decomposition, amplitude equations and supported divergence-free oscillatory waves."},
    {"id": 5, "name": "compact mean corrections", "status": "pending",
     "remaining": "Section 8 compact mean corrections and the associated defect solve."},
    {"id": 6, "name": "residual-improvement iteration", "status": "pending",
     "remaining": "Section 9 recursive choices, residual-improvement cycle, summation and convergence control."},
    {"id": 7, "name": "final compact field and force", "status": "formal-structure",
     "implemented": (
         "Section 10 support/plateau geometry, explicit C-infinity cutoff representative, analytic cutoff gradient wired into "
         "curl(cA)=c curl(A)+grad(c) cross A, and numerical forcing reconstruction/independent manufactured diagnostics."
     ),
     "remaining": (
         "Justify the executable transition-collar representative, implement paper-specific time activation and smooth extension through t=1, "
         "and connect completed upstream local-field inputs."
     )},
    {"id": 8, "name": "independent verification", "status": "diagnostic-only",
     "implemented": "Manufactured solutions, refinement checks, divergence/energy diagnostics, labelled exports and fail-closed provenance audit.",
     "remaining": "Validate an actual paper instance and add uniform/interval or formal full-field certificates."},
]


def construction_status() -> dict:
    return {
        "schema_version": 1,
        "paper_exact_velocity_available": False,
        "status": "partial-executable-reconstruction",
        "remote_publication_performed": True,
        "limitations": [
            "The bundled Gaussian profile is a toy, not OpenAI's constructed profile.",
            "Numerical residual tests are not the paper's proof or a Lean certificate.",
            "A force reconstructed from the same residual stencil is not independent evidence.",
            "Pointwise or sampled margins are not uniform-in-parameter certificates.",
            "The executable Section 10 transition collar has the official support/plateau geometry but is not claimed pointwise identical to Mathlib's noncomputable ContDiffBump.",
        ],
        "sources": copy.deepcopy(SOURCE_PINS),
        "stages": copy.deepcopy(STAGES),
    }
