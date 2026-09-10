"""Machine-readable scope and source pins; test success cannot change truth status."""
from __future__ import annotations
import copy

SOURCE_PINS = {
    "repository": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction",
    "baseline_commit": "b0d964ab4a69c25863078830faf76296550ae4ea",
    "integration_base_commit": "b01aaeebc6ec8a39f6692109b6a3c81b92cd0f32",
    "integration_base_tree_verified": "8788a0faf709c02ece63a41774d29f9754d7a59b",
    "paper_url": "https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf",
    "paper_equations_checked": ["4.1-4.7", "5.1", "5.27", "Lemma A.2"],
    "upstream_lean_repository": "openai/NavierStokesAndEuler",
    "upstream_commit_metadata_observed": "f9e8bc5b38b6e212696e8a30e3e91517af887bbd",
    "lean_source_reviewed": False,
    "lean_compiled": False,
    "paper_hash_verified": False,
}
STAGES = [
    {"id":0,"name":"similarity geometry","status":"implemented-numerically",
     "remaining":"No interval or arbitrary-precision certificate."},
    {"id":1,"name":"leading profile","status":"partial",
     "implemented":"Kinematics, pressure interface, pointwise moment-solver primitives.",
     "remaining":"Construct Theorem 4.6 profiles, matching data, cone and moment constraints."},
    {"id":2,"name":"all-order background","status":"formal-structure",
     "implemented":"Finite analytic cutoff/curl assembly from supplied coefficients.",
     "remaining":"Coefficient recursion, admissible cutoff schedule, all-order estimates."},
    {"id":3,"name":"dyadic charts and transported phases","status":"pending",
     "remaining":"Implement Section 6 choices and transported wave data."},
    {"id":4,"name":"oscillatory stress realization","status":"pending",
     "remaining":"Implement Section 7 stress cone and amplitudes."},
    {"id":5,"name":"compact mean corrections","status":"pending",
     "remaining":"Section 8 radial inverses and five-equation defect solve."},
    {"id":6,"name":"residual-improvement iteration","status":"pending",
     "remaining":"Section 9 recursive choices, summation and convergence control."},
    {"id":7,"name":"final compact field and force","status":"formal-structure",
     "implemented":"Localization product rule and independent residual diagnostics.",
     "remaining":"Instantiate the field and force smooth through t=1; no force extension supplied."},
    {"id":8,"name":"independent verification","status":"partial",
     "implemented":"Regression, manufactured solutions, refinement, finite-cylinder energy, labelled data export.",
     "remaining":"Validate an actual paper instance, uniform estimates and formal cross-checks."},
]


def construction_status() -> dict:
    return {
        "schema_version":1,
        "paper_exact_velocity_available":False,
        "status":"partial-executable-reconstruction",
        "remote_publication_performed":False,
        "limitations":[
            "The bundled Gaussian profile is a toy, not OpenAI's constructed profile.",
            "Numerical residual tests are not the paper's proof or a Lean certificate.",
            "R(u,p) minus the same computed R(u,p) is not an independent forcing check.",
            "Pointwise moment smallness is not a uniform-in-parameter bound.",
        ],
        "sources":copy.deepcopy(SOURCE_PINS),
        "stages":copy.deepcopy(STAGES),
    }
