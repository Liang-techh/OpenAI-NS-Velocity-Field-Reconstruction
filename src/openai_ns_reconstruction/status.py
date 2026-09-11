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
        "5.1-5.7",
        "5.10-5.16",
        "5.27",
        "6.1-6.11",
        "7.1-7.11",
        "10.4",
        "Section 10 time localization",
        "Section 10 closed-past extension",
        "Section 10 endpoint Taylor-Borel extension",
        "Lemma 5.1",
        "Lemma 5.2",
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
         "NaturalAxisRange h/j bounds, unique H-root bracket, delta=j/10, an explicit B=2 ideal-prefix "
         "PressureDatum.Admissible witness with closed-form axis pressure, the pinned Gaussian-flat outgoing scalar "
         "schedule with analytic sigma', an explicit S=32 derivative-bound witness, constructive flattenLength, "
         "SchedulePressure.shapeExponent, the complete executable OutgoingTail finalAngular(y,eta) / "
         "SchedulePressure.clockWeight(y) chain, the actual all-real-line SchedulePressure.axisPressure evaluator, "
         "an analytic no-sampling low-|Z| uniform H^2 margin that instantiates the theorem-side sigma=sqrt(m)/20 choice, "
         "and the theorem-faithful Lambda/C selector Lambda=max(1+B+L,1+(14000/9)B), "
         "C=exp(Lambda*realPartSup), conditional on certified coefficient-family B/L and complex phase realPartSup."
     ),
     "remaining": (
         "Materialize the actual coefficient-family remainderBound/remainderLip and certify the complex compact-set "
         "realPartSup used by the landed Lambda/C selector; use those certified constants to materialize the "
         "coefficient-space fixed-point phi/u/average/pressure fields, connect them to NaturalProfileAssembly, then "
         "verify the paper's support, moments, matching and cone conditions."
     )},
    {"id": 2, "name": "all-order background", "status": "formal-structure",
     "implemented": (
         "Finite coefficient assembly, Eq. (5.2) radial flux, Eq. (5.27) streamfunction/vector potential, analytic "
         "cutoff/curl terms, the Eq. (5.5) pressure recurrence row with Pi_n(0,eta)=0 integration, the regular Eq. (5.6) "
         "Omega_k/X source evaluator from second jets of V_j=X v_j without axis division, the Eq. (5.7) six-component "
         "singular inverse G with one-step Picard map plus the first positive-order G f_n for supplied A0/A1/f_n, the "
         "Lemma 5.2 / Eqs. (5.14)-(5.16) compact five-moment repair using two U bumps and three E bumps, the "
         "Eq. (5.15) forward reconstruction of F_n/V_n/Pi_n from supplied repaired analytic coefficient data using "
         "the axis-regular Eq. (5.2)/(5.5) paths, and a finite-prefix SlowBorelBase/DiagonalScale recursive cutoff-scale "
         "constructor from supplied analytic normalized-template bounds C[j,m], with exact edge certificates and the doubling envelope."
     ),
     "remaining": (
         "Derive the paper's profile-dependent A0/A1/f_n from materialized leading/lower-order data, iterate Eq. (5.7) "
         "to the Lemma 5.1 solution and connect it to the constructed Eq. (5.6) source; materialize the true eta-dependent "
         "repaired coefficient hierarchy and certify the support/stress conclusions required before reuse at the next order; "
         "derive the true uniform compactness bounds C[j,m] from that hierarchy and instantiate the full recursive cutoff "
         "sequence; then prove local finiteness and arbitrary-order truncation/residual decay."
     )},
    {"id": 3, "name": "dyadic charts and transported phases", "status": "formal-structure",
     "implemented": (
         "Fixed dyadic geometry, an Eq. (6.8) active-shell/slow-label bridge with q/Q and X reconstruction, S_*^-3 "
         "slow-support enclosure and typed TangentialBaseJetProvider freezing Eq. (6.11)/(7.2) representative data, "
         "paper-admissible normalized C-infinity translate squared partitions for dyadic q and the S_*^-3 (R,Z,T) "
         "product mesh with sign-duplicated labels sharing one cutoff, a constructive Section 6.2 auxiliary-slot witness "
         "with the pinned 2250-color mod-9/mod-5/sign data, exact rational centers, covering-index checks and an explicit "
         "common conservative r0 under the discrete interaction hypotheses, a physical slow-support adjacency bridge "
         "using the axis-dependent Q^(1/2), Q^(1/2-h), Q widths and a common physical-point cross-band certificate, "
         "Section 7.1 phase/wavevector/tangent-frame algebra, the pinned PhaseEstimates scalar scale gate for the "
         "Eq. (7.9)-(7.11) error envelope, and a fail-closed pointwise rounded_normal_estimates adapter using the pinned "
         "Lean roundedFrequency semantics."
     ),
     "remaining": (
         "Instantiate the TangentialBaseJetProvider from the paper-exact Proposition 5.5 background and theorem choices, "
         "and certify uniform LocalBaseBounds/C1-C2 hypotheses on every active box so the pointwise adapter becomes a "
         "uniform phase/frame/damping certificate; then continue to the stress-cone, amplitude/curl wave, mean-correction "
         "and residual-improvement stages."
     )},
    {"id": 4, "name": "oscillatory stress realization", "status": "pending",
     "remaining": "Stress-cone decomposition, amplitude equations and supported divergence-free oscillatory waves."},
    {"id": 5, "name": "compact mean corrections", "status": "pending",
     "remaining": "Section 8 compact mean corrections and the associated defect solve."},
    {"id": 6, "name": "residual-improvement iteration", "status": "pending",
     "remaining": "Section 9 recursive choices, residual-improvement cycle, summation and convergence control."},
    {"id": 7, "name": "final compact field and force", "status": "formal-structure",
     "implemented": (
         "Section 10 support/plateau geometry, explicit C-infinity spatial cutoff representative, analytic cutoff gradient "
         "wired into curl(cA)=c curl(A)+grad(c) cross A, pinned time-switch support/plateau geometry with analytic chi', "
         "activated velocity/pressure adapters and an independent numerical cross-check of the TimeLocalization residual "
         "identity, the pinned closed-past zeroBefore/pastVelocity/pastPressure branch plus a diagnostic evaluation of the "
         "actual closed-past NS residual on 0<t<1, the pinned endpoint Taylor-Borel algebra with exact "
         "DiagonalScale.doublingEnvelope recurrence and locally finite right-extension evaluation, a fail-closed dense "
         "full-spacetime endpoint-jet adapter that contracts CandidateFromLimits tensors in every derivative slot with "
         "timeVector=(1,0), and an exact-rational SpatialBorelExtension boundSum/localScale/doubling scale-schedule adapter "
         "with supplied-bound 2^-j derivative-tail certificates."
     ),
     "remaining": (
         "Justify the executable spatial/time transition-collar representatives where the formalization uses noncomputable "
         "ContDiffBump values; materialize the actual closed-past localized NS residual's full spacetime derivative family "
         "and prove all locally uniform t->1- limits; derive genuine analytic compact-template derivative bounds for those "
         "true jets and instantiate the landed scale schedule; then prove jet matching and smooth global force extension "
         "through t=1; also connect completed upstream local-field inputs."
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
            "The executable Section 10 spatial/time transition collars have the official support/plateau geometry but are not claimed pointwise identical to Mathlib's noncomputable ContDiffBump.",
            "The endpoint Borel/full-spacetime-jet/scale-schedule adapters still depend on caller-supplied endpoint jets or analytic template bounds and do not certify the actual residual derivative limits or smooth t=1 force extension.",
        ],
        "sources": copy.deepcopy(SOURCE_PINS),
        "stages": copy.deepcopy(STAGES),
    }
