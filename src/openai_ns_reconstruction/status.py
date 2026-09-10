"""Machine-readable scope and source pins; passing tests cannot change truth status."""
from __future__ import annotations
import copy

SOURCE_PINS = {
    "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction",
    "integration_base_commit": "e832242c6ab389b602037659e719cd2b87cb7d2c",
    "integrated_pr5_head": "9c7ad1ff6dd4b21806ff5ee2248974975542b76d",
    "paper_url": "https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf",
    "paper_equations_checked": ["4.1-4.7", "5.1-5.2", "5.27", "6.1-6.6", "7.1-7.8", "10.4", "Lemma A.2"],
    "upstream_lean_repository": "openai/NavierStokesAndEuler",
    "upstream_commit_metadata_observed": "f9e8bc5b38b6e212696e8a30e3e91517af887bbd",
    "lean_source_reviewed": False,
    "lean_compiled": False,
    "paper_hash_verified": False,
}

STAGES = [
    {"id": 0, "name": "similarity geometry", "status": "implemented-numerically",
     "implemented": "Stable direct-tau chart plus analytic coordinate derivatives.",
     "remaining": "No interval or arbitrary-precision certificate."},
    {"id": 1, "name": "leading profile", "status": "partial",
     "implemented": "Kinematics, pressure interface, natural-axis rescaling, heat exterior, pointwise moment primitives.",
     "remaining": "Materialize the complete Theorem 4.6 fixed-point profile and existential parameter choices."},
    {"id": 2, "name": "all-order background", "status": "formal-structure",
     "implemented": "Finite coefficient assembly, Eq. (5.2) radial flux, analytic cutoff/curl terms.",
     "remaining": "Recursive coefficient solver, moment repairs, certified cutoff schedule and all-order estimates."},
    {"id": 3, "name": "dyadic charts and transported phases", "status": "formal-structure",
     "implemented": "Dyadic chart geometry and Section 7.1 phase/frame formulas for supplied base jets.",
     "remaining": "Instantiate paper-exact base fields and verify uniform phase/frame/damping estimates."},
    {"id": 4, "name": "oscillatory stress realization", "status": "pending",
     "remaining": "Stress-cone amplitudes and supported divergence-free oscillatory waves."},
    {"id": 5, "name": "compact mean corrections", "status": "pending",
     "remaining": "Section 8 radial inverses and five-equation defect solve."},
    {"id": 6, "name": "residual-improvement iteration", "status": "pending",
     "remaining": "Section 9 recursive choices, summation and convergence control."},
    {"id": 7, "name": "final compact field and force", "status": "formal-structure",
     "implemented": "Section 10 support geometry, localization product rule and numerical forcing reconstruction.",
     "remaining": "Completed local field inputs and smooth forcing extension through t=1."},
    {"id": 8, "name": "independent verification", "status": "partial",
     "implemented": "Manufactured solutions, Taylor-Green refinement, divergence, finite-cylinder energy, labelled exports.",
     "remaining": "Validate an actual paper instance, uniform estimates and formal/interval cross-checks."},
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
            "Pointwise moment smallness is not a uniform-in-parameter bound.",
        ],
        "sources": copy.deepcopy(SOURCE_PINS),
        "stages": copy.deepcopy(STAGES),
    }
