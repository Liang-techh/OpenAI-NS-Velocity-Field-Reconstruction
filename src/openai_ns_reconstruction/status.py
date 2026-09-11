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
        "5.1-5.8",
        "5.10-5.16",
        "5.27",
        "6.1-6.11",
        "7.1-7.11",
        "7.24-7.30",
        "10.4",
        "10.20-10.21",
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
         "the theorem-faithful Lambda/C selector Lambda=max(1+B+L,1+(14000/9)B), C=exp(Lambda*realPartSup), a "
         "fail-closed AxisContraction-style bound/Lipschitz propagation from certified coefficient/operator norm ledgers "
         "to conservative remainderBound/remainderLip values, the canonical analytic-input norm adapter: from a "
         "certified common complex neighborhood radius rho>0 and common field-value bound B it uses epsilon=rho/2, "
         "radiusLoss(1/2)=12, fixed-field norm <=12B, normalized amplitude M=12, ||chi||<=12B, K<=30720B, and a "
         "conservative complete factorial-series enclosure for the naturalResolvent norm, plus an actual-schedule "
         "analytic-neighborhood certificate that derives a positive common rho, the eleven-field common complex bound B, "
         "and a conservative axisPhase real-part supremum from theorem-side sigma and analytic inequalities without "
         "grid maxima or fitted complex data. The actual SchedulePressure certificate is now wired through the "
         "coefficient-family norm / AxisResolvent / remainder / Lambda-C chain without caller-supplied rho, B, K, "
         "resolvent norm or scale constants; for the current theorem-admissible regression schedule a downward-rounded "
         "positive factorial-majorant term already exceeds binary64 range, so stage1_scale_chain fails closed before "
         "inventing remainder or Lambda/C values."
     ),
     "remaining": (
         "Tighten the connected actual-schedule rho/B/realPartSup certificate-to-contraction chain, or carry its "
         "factorial majorant in a theorem-faithful wider representation, while preserving every contraction hypothesis; "
         "the current conservative regression chain is blocked by a certified binary64 representation overflow, not by "
         "a lower bound on the true resolvent. Then materialize the coefficient-space fixed-point phi/u/average/pressure "
         "fields, connect those genuine fields to NaturalProfileAssembly, and verify the paper's support, moments, "
         "matching and cone conditions. The current binary64 analytic inequalities are executable conservative "
         "certificates, not interval or Lean proof objects."
     )},
    {"id": 2, "name": "all-order background", "status": "formal-structure",
     "implemented": (
         "Finite coefficient assembly, Eq. (5.2) radial flux, Eq. (5.27) streamfunction/vector potential, analytic "
         "cutoff/curl terms, the Eq. (5.5) pressure recurrence row with Pi_n(0,eta)=0 integration, the regular Eq. (5.6) "
         "Omega_k/X source evaluator from second jets of V_j=X v_j without axis division, the Eq. (5.7) six-component "
         "singular inverse G with one-step Picard map plus the first positive-order G f_n for supplied A0/A1/f_n, a "
         "theorem-shaped Lemma 5.1 / Eq. (5.8) Picard-term majorant and complete-tail truncation certificate using the "
         "p_k=ceil(k/2) Cauchy-loss factor and a fail-closed two-step geometric tail gate, the Lemma 5.2 / "
         "Eqs. (5.14)-(5.16) compact five-moment repair using two U bumps and three E bumps, the Eq. (5.15) forward "
         "reconstruction of F_n/V_n/Pi_n from supplied repaired analytic coefficient data using the axis-regular "
         "Eq. (5.2)/(5.5) paths, a finite-prefix SlowBorelBase/DiagonalScale recursive cutoff-scale constructor from "
         "supplied analytic normalized-template bounds C[j,m] with exact edge certificates and the doubling envelope, an "
         "exact finite-prefix cutoff support certificate recording plateau, the unique possible transition order, "
         "forced-zero tail within the constructed prefix and stable truncation once that zero tail is reached, a "
         "fixed-prefix tail-order arithmetic gate recording the Lean bounds 2^-J q^(h(J+1)-m) and physical exponent "
         "h(J+1)+b-2M, solving exactly for the minimal J needed for a requested finite decay order and failing closed when "
         "the constructed finite scale schedule is too short, and an exact SlowExpansionResidual finite recurrence / "
         "truncation bridge that keeps every nonzero retained recurrence defect explicit and, only after independent exact "
         "recurrence cancellation, bounds the omitted pair/shift tail by the first omitted slow factor "
         "q^(b+2(N+1)h) times its coefficient majorant."
     ),
     "remaining": (
         "Derive the paper's profile-dependent A0/A1/f_n from materialized leading/lower-order data, verify the actual "
         "analytic constants entering the landed Eq. (5.8) convergence/tail certificate, and connect the converged "
         "Eq. (5.7) Picard series to the Eq. (5.6) source and solved coefficient hierarchy. Materialize the true "
         "eta-dependent repaired hierarchy and certify the support/stress conclusions required before reuse at the next "
         "order; derive true uniform C[j,m] bounds, instantiate the full infinite recursive cutoff sequence, then prove "
         "theorem-level local finiteness and the actual coefficient identities cancelling the NS residual order by order. "
         "Only those independently proved identities may be fed to the landed finite recurrence/truncation bridge before "
         "claiming Proposition 5.3/all-jets-flat residual decay."
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
         "Eq. (7.9)-(7.11) error envelope, a fail-closed pointwise rounded_normal_estimates adapter using the pinned Lean "
         "roundedFrequency semantics, a fail-closed UniformLocalBase bridge from independently certified LocalBase sup "
         "bounds plus slow-box diameter d<=S^-3 to the exact M(S^-3+epsilon^2) derivative-error envelope, and a "
         "theorem-shaped BasePhaseGeometry frame/damping implication layer: it preserves the distinction between "
         "normalConstant(M) and the enlarged family phaseConstant(M)=normalConstant(frequencyBound(M)), bounds the three "
         "moving-frame coefficient errors by 16 G^2(1+G)E and the viscosity discrepancy by 4 M(2A+5)delta, and fails "
         "closed on its scalar large-band hypotheses without pretending the missing vector/analytic assumptions are proved."
     ),
     "remaining": (
         "Instantiate the TangentialBaseJetProvider from the paper-exact Proposition 5.5 background and theorem choices, "
         "and analytically certify the actual LocalBaseBounds/C1-C2 hypotheses on every active box together with the unit, "
         "orthogonality, normal-closeness and slot-normal derivative-closeness assumptions required by BasePhaseGeometry; "
         "then feed those genuine hypotheses through the landed uniform normal and frame/damping envelopes before using the "
         "stress layer on actual pulse data."
     )},
    {"id": 4, "name": "oscillatory stress realization", "status": "formal-structure",
     "implemented": (
         "Pinned Section 7 signed two-slot reference covariance algebra for the target (-m,t): exact determinant, strict "
         "positive-cone condition |a t|<b m, explicit positive squared amplitudes, the manuscript ratio specialization "
         "a=-c_*sqrt(1+u_*^2) with |c_* t/m|<u_*/sqrt(1+u_*^2), and the exact epsilon*mask^2 covariance scaling after "
         "sqrt(epsilon)*mask amplitude localization. A separate fail-closed actual-vs-reference 2x2 covariance "
         "perturbation certificate now accepts independently certified normalized column errors, derives a positive "
         "determinant margin, a conservative spectral ||H^-1||_2 upper bound and a positive lower bound for the perturbed "
         "squared amplitudes, with independent NumPy determinant/inverse/solve cross-checks on synthetic fixtures."
     ),
     "remaining": (
         "Construct the actual Eq. (7.27) pulse-integrated covariance columns and prove the uniform Eq. (7.28) analytic "
         "column-error bound from the real background/phase/damping data; feed those genuine errors through the landed "
         "determinant/inverse stability and positive-amplitude certificate for the real T_{0,*} data, then solve the "
         "amplitude transport/ODE and realize the oscillatory corrections by supported divergence-free curls."
     )},
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
         "timeVector=(1,0), an exact-rational SpatialBorelExtension boundSum/localScale/doubling scale-schedule adapter "
         "with supplied-bound 2^-j derivative-tail certificates and an exact omitted derivative-series certificate "
         "sum_{j>N} 2^-j = 2^-N whenever N>=max(spatial window, derivative order), a value-level CandidateFromLimits "
         "trace/glue bridge that preserves the supplied closed-past residual for t<T, uses only the validated degree-zero "
         "endpoint tensor for t>=T, and composes the normal-jet contraction with the existing locally finite Borel future "
         "branch, an exact support-cylinder implication from SpatialLocalization: support lies in r^2<=1/16, |z|<=1/4 "
         "(volume pi/32), so an independently certified fixed-time |u|<=M bound yields E(t)<=pi M^2/64 without inferring "
         "M from samples, and a fail-closed late-origin preservation certificate for 3/4<=t<1 showing the official "
         "spatial plateau gives c=1 and grad(c)=0 at the origin while the time switch is exactly one, so final localization "
         "preserves any upstream-certified origin curl value without sampling a blow-up sequence."
     ),
     "remaining": (
         "Justify the executable spatial/time transition-collar representatives where the formalization uses noncomputable "
         "ContDiffBump values; materialize the actual closed-past localized NS residual's full spacetime derivative family "
         "and prove all locally uniform t->1- limits; derive genuine analytic compact-template derivative bounds for those "
         "true jets and instantiate the landed scale schedule and exact 2^-N derivative-tail budget; then prove the traced "
         "degree-zero endpoint value and all normal jets arise from those actual limits so the resulting global force is "
         "smooth through t=1. Independently establish the paper's uniform bounded kinetic-energy conclusion along the "
         "actual blow-up limit, connect completed upstream local-field inputs, and prove the required upstream "
         "||curl A(t,0)||->infinity premise; the landed late-origin adapter preserves that premise but does not establish it."
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
            "The endpoint Borel/full-spacetime-jet/scale-schedule/traced-residual adapters still depend on caller-supplied endpoint jets or analytic template bounds and do not certify the actual residual derivative limits or smooth t=1 force extension; the exact 2^-N omitted-derivative tail budget is conditional on those genuine bounds.",
            "The Section 10 compact-support energy adapter is only a fixed-time implication from an independently certified sup bound; it does not prove a uniform bounded-energy estimate as t approaches the blow-up time.",
            "The Section 10 late-origin adapter preserves an upstream-certified origin curl value on 3/4<=t<1 but does not prove that the upstream curl actually blows up.",
            "The Section 7 stress-cone module is only the exact reference algebra; the separate covariance-perturbation module is only a caller-certified finite-dimensional implication, while actual pulse integrals, the Eq. (7.28) analytic error estimate, amplitude transport and wave realization are still missing.",
            "The connected Stage-1 actual-schedule scale chain currently stops on a certified binary64 overflow of a conservative positive resolvent-majorant term; this is not a lower bound on the true resolvent or a failure of the existence theorem.",
        ],
        "sources": copy.deepcopy(SOURCE_PINS),
        "stages": copy.deepcopy(STAGES),
    }
