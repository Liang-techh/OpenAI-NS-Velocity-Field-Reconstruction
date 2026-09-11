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
         "resolvent norm or scale constants. A theorem-faithful componentwise refinement applies the same "
         "radiusLoss(1/2)=12 estimate to each of the eleven field bounds B_k, keeps max_k 12B_k only as the "
         "CoefficientFamily.bound, and uses the chi-specific ||chi||<=12B_chi in AxisResolvent. On the current "
         "theorem-admissible regression schedule this removes the earlier cross-field resolvent overflow: the complete "
         "factorial resolvent majorant is representable in binary64. The original binary64 remainder path still fails "
         "closed during conservative remainderBound/remainderLip propagation, but the landed wide-arithmetic adapter now "
         "propagates the same positive Controlled/remainder majorant algebra with 96-digit Decimal arithmetic rounded "
         "toward +infinity. On that actual regression datum its remainderBound/remainderLip ledger remains finite even "
         "when a conservative bound exceeds sys.float_info.max. The landed wide Lambda/C selector now consumes those "
         "Decimal bounds without binary64 down-conversion, evaluates Lambda=max(1+B+L,1+(14000/9)B) with upward-rounded "
         "Decimal arithmetic, and represents C=exp(Lambda*phaseSup) by its certified symbolic exponent so an enormous but "
         "finite theorem threshold is not mistaken for numerical infinity. The actual-schedule stage1_scale_chain_wide "
         "path supplies no caller-chosen rho, field bounds, resolvent, remainder, Lambda, or C. The landed "
         "axis_fixed_point_picard adapter then carries this actual-schedule wide chain through the pinned scalar "
         "AxisContraction gates s*remainderBound<=1 and s*remainderLip<=1/2 with s=1/(2Lambda), using upward-rounded "
         "Decimal arithmetic, and exposes the theorem-shaped Picard update x_(n+1)=x0+s*naturalRemainder(x_n) together "
         "with a geometric tail budget. It deliberately leaves coefficient-state operations and the genuine "
         "naturalRemainder evaluator as typed dependencies; therefore the gate is executable but no coefficient-space "
         "fixed point has been materialized."
     ),
     "remaining": (
         "Implement a genuine compatible AxisCoefficientSpace/coefficientOperators/naturalRemainder backend, execute the "
         "landed theorem-shaped Picard map from the pinned referencePair on actual coefficient states, and materialize the "
         "coefficient-space fixed-point phi/u fields from the actual schedule/operator/remainder data. Derive the "
         "corresponding average/pressure fields and connect those genuine fields together with the landed wide Lambda/C "
         "choice to NaturalProfileAssembly. If a practical executable solve needs smaller constants, tighten the theorem-"
         "side conservative majorants without weakening the contraction hypotheses; the former binary64 Lambda/C "
         "representation blocker is closed, and the historical overflow is not a lower bound on the true "
         "resolvent/remainder or a failure of the Lean existence theorem. Then verify the paper's support, moments, "
         "matching and cone conditions. The current analytic inequalities, wide arithmetic and contraction gate are "
         "executable conservative certificates, not interval or Lean proof objects."
     )},
    {"id": 2, "name": "all-order background", "status": "formal-structure",
     "implemented": (
         "Finite coefficient assembly, Eq. (5.2) radial flux, Eq. (5.27) streamfunction/vector potential, analytic "
         "cutoff/curl terms, the Eq. (5.5) pressure recurrence row with Pi_n(0,eta)=0 integration, the regular Eq. (5.6) "
         "Omega_k/X source evaluator from second jets of V_j=X v_j without axis division, the Eq. (5.7) six-component "
         "singular inverse G with one-step Picard map plus the first positive-order G f_n for supplied A0/A1/f_n. A "
         "landed PositiveAxisSystem adapter now constructs the pinned real A0, sparse A1 and f_n at lambda_n=2nh from "
         "typed BaseJet(phi,axial,beta) and SourceJet(angular,axial,pressureProduct,omegaQuotient) providers, so the first "
         "positive-order G f_n can be driven by theorem-shaped coefficient data rather than an arbitrary six-vector. A "
         "theorem-shaped Lemma 5.1 / Eq. (5.8) Picard-term majorant and complete-tail truncation certificate uses the "
         "p_k=ceil(k/2) Cauchy-loss factor and a fail-closed two-step geometric tail gate, the Lemma 5.2 / "
         "Eqs. (5.14)-(5.16) compact five-moment repair using two U bumps and three E bumps, and a finite eta-jet "
         "adapter that propagates supplied unrepaired moment derivatives and p=e_*f derivatives through exact "
         "Leibniz/quotient recurrences and the same fixed B_U/B_E repair matrices without pointwise fitting. A landed "
         "function-level compact-repair adapter now turns those coefficient jets into actual eta-dependent repaired "
         "U_n(X,eta), E_n(X,eta), and d_eta U_n functions, preserving the base profile outside the five bump supports; "
         "its LeadingProfile bridge hard-codes paper_exact=False and never silently inherits an unrepaired pressure. The "
         "Eq. (5.15) forward reconstruction of F_n/V_n/Pi_n from supplied repaired analytic coefficient data uses the "
         "axis-regular Eq. (5.2)/(5.5) paths. Also present are a finite-prefix SlowBorelBase/DiagonalScale recursive "
         "cutoff-scale constructor from supplied analytic normalized-template bounds C[j,m] with exact edge certificates "
         "and the doubling envelope. When a theorem-admissible local scale exceeds binary64, the constructor keeps the "
         "same dyadic inequality and produces an arbitrary-precision power-of-two Python integer witness, checking the "
         "original inequality in log space rather than clipping the scale; reciprocal_support_log_edge preserves the "
         "support edge in log form, while the binary64 reciprocal accessor fails closed if it would underflow to fake zero. "
         "An exact finite-prefix cutoff support certificate records plateau, the unique possible transition order, forced-"
         "zero tail within the constructed prefix and stable truncation once that zero tail is reached, a fixed-prefix "
         "tail-order arithmetic gate records the Lean bounds 2^-J q^(h(J+1)-m) and physical exponent h(J+1)+b-2M, "
         "solving exactly for the minimal J needed for a requested finite decay order and failing closed when the "
         "constructed finite scale schedule is too short, and an exact SlowExpansionResidual finite recurrence / "
         "truncation bridge keeps every nonzero retained recurrence defect explicit and, only after independent exact "
         "recurrence cancellation, bounds the omitted pair/shift tail by the first omitted slow factor "
         "q^(b+2(N+1)h) times its coefficient majorant."
     ),
     "remaining": (
         "Instantiate the paper-derived A0/A1/f_n adapter from materialized leading and lower-order BaseJet/SourceJet "
         "data, verify the actual C_n, analytic radii and Delta entering the landed Eq. (5.8) convergence/tail certificate, "
         "and connect the converged Eq. (5.7) Picard series to the Eq. (5.6) source and solved coefficient hierarchy. "
         "Materialize the true eta-dependent repaired hierarchy by deriving the unrepaired moment functions, genuine "
         "moment/p=e_*f jets and a certified nonvanishing patch factor from that hierarchy, then certify support/stress "
         "conclusions before reuse at the next order; derive true uniform C[j,m] bounds, instantiate the full infinite "
         "recursive cutoff sequence, then prove theorem-level local finiteness and the actual coefficient identities "
         "cancelling the NS residual order by order. Only those independently proved identities may be fed to the landed "
         "finite recurrence/truncation bridge before claiming Proposition 5.3/all-jets-flat residual decay."
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
         "perturbation certificate accepts independently certified normalized column errors, derives a positive "
         "determinant margin, a conservative spectral ||H^-1||_2 upper bound and a positive lower bound for the perturbed "
         "squared amplitudes, with independent NumPy determinant/inverse/solve cross-checks on synthetic fixtures. A "
         "theorem-shaped pulse-covariance error-budget adapter now turns independently certified Lemma-7.4 ratio/frame "
         "error E_ratio and positive-weight normalized first-moment bound M1 into "
         "|e_sigma|<=E_ratio+u_*sqrt(1+c_*^2)M1, and collapses certified O(S_*^-1)+O(S_*^-1/2) inputs to the Eq. (7.28) "
         "C/sqrt(S_*) envelope without identifying synthetic/caller data with the actual pulse integral. The coefficient-"
         "level curl-realization algebra from the pinned Formula (30) is also executable: B=|n|^-2(n cross a), the "
         "principal harmonic is the tangential projection, inverseCarrier=i/K, and the displayed derivative remainder is "
         "(i/K)curl(B). Exact tangency n.a=0 is required before exposing the tangent specialization; caller-supplied "
         "coefficient_curl is never treated as proof of the genuine cylindrical curl. The landed theorem-shaped primary "
         "amplitude ODE algebra now maps MovingFrameODE/GrowingMode/PrimaryODE pointwise: moving-frame a/b/c, reference "
         "defects B=b-lambda/h and C=c-lambda*h, four modal errors, j^2 viscosity damping, the exact 2x2 modal operator, "
         "projected forcing transform, and x=p+q, y=h(p-q) with h'=rate*h. Its regression independently evaluates the "
         "physical rhsX/rhsY path and differentiates the modal basis. A separate finite-interval PrimaryODE adapter now "
         "integrates the supplied time-dependent modal IVP, reconstructs the moving physical amplitudes, and independently "
         "checks the Volterra integral defect with quadrature; regression uses a closed-form constant-coefficient oracle, "
         "not the production solver itself. This remains formal-structure execution for supplied datum/forcing and is not "
         "the paper's Proposition 5.5/pulse instance or a rigorous numerical error certificate."
     ),
     "remaining": (
         "Construct the actual pulse-integrated covariance columns of Eq. (7.27), including the true pulse functions, "
         "from the paper-exact background; prove the genuine Lemma 7.4 ratio/frame bound and Gaussian normalized "
         "first-moment estimate on every active pulse, then feed those certified inputs through the landed Eq. (7.28) "
         "error-budget and perturbation certificates to obtain actual-data determinant/inverse stability and positive-cone "
         "persistence for the real T_{0,*} data. Instantiate the landed primary-amplitude datum and forcing from the genuine "
         "Proposition 5.5/phase-frame/pulse fields, use the landed finite-interval PrimaryODE Volterra solution adapter on "
         "those actual inputs, and establish its ODE/initial-value/weighted estimates with an appropriate rigorous error "
         "certificate. Only then connect that actual amplitude to the landed curl algebra with the actual localized "
         "coefficient, cylindrical derivative, smoothness, support, zero-germ and phase-patch hypotheses required for a "
         "genuine supported divergence-free oscillatory wave."
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
         "DiagonalScale.doublingEnvelope recurrence and locally finite right-extension evaluation, including a fail-closed "
         "close_left_open_past adapter for upstream evaluators defined only on t<T: it leaves t<T unchanged, fills t=T "
         "from jet(0,x), rejects t>T, and glue_to_left_open_past preserves the pinned closed-past branch at the join without "
         "asserting continuity. A fail-closed dense full-spacetime endpoint-jet adapter contracts CandidateFromLimits "
         "tensors in every derivative slot with timeVector=(1,0), an exact-rational SpatialBorelExtension "
         "boundSum/localScale/doubling scale-schedule adapter has supplied-bound 2^-j derivative-tail certificates and an "
         "exact omitted derivative-series certificate sum_{j>N} 2^-j = 2^-N whenever N>=max(spatial window, derivative "
         "order), a value-level CandidateFromLimits trace/glue bridge preserves the supplied closed-past residual for t<T, "
         "uses only the validated degree-zero endpoint tensor for t>=T, and composes the normal-jet contraction with the "
         "existing locally finite Borel future branch, an exact support-cylinder implication from SpatialLocalization "
         "gives support in r^2<=1/16, |z|<=1/4 (volume pi/32), so an independently certified fixed-time |u|<=M bound "
         "yields E(t)<=pi M^2/64 without inferring M from samples, and a fail-closed late-origin preservation certificate "
         "for 3/4<=t<1 shows the official spatial plateau gives c=1 and grad(c)=0 at the origin while the time switch is "
         "exactly one, so final localization preserves any upstream-certified origin curl value without sampling a blow-up "
         "sequence. A separate fail-closed endpoint-limit majorant records the elementary sufficient condition "
         "sup_x||partial_t D^n R||<=C(T-t)^(-alpha), 0<=alpha<1, and derives the locally-uniform Cauchy/endpoint-tail "
         "modulus C/(1-alpha)(T-t)^(1-alpha) without fitting C or alpha from samples. A landed Section 10 endpoint-"
         "localization transfer now accepts this majorant only on the official late time plateau: it requires T=1 and a "
         "majorant start time at or after 3/4, then checks timeSwitch=1, timeSwitch'=0 and chi^2-chi=0 there before "
         "transferring the upstream residual derivative bound unchanged toward the endpoint. Transition-collar or non-T=1 "
         "majorants fail closed. A landed pointwise endpoint-force support certificate follows "
         "SpatialBorelExtension.extension_zero_of_coefficients_zero outside the official cylinder: it checks the past "
         "residual exactly for t<T, the degree-zero endpoint trace at t=T, and exactly the locally finite active future "
         "normal coefficients for t>T, requiring exact zeros with no tolerance or sampled support inference."
     ),
     "remaining": (
         "Justify the executable spatial/time transition-collar representatives where the formalization uses noncomputable "
         "ContDiffBump values; materialize the actual closed-past localized NS residual's full spacetime derivative family "
         "and prove genuine derivative bounds (or equivalent estimates) on the official t>=3/4 plateau that discharge the "
         "landed endpoint-limit majorant through the fail-closed T=1 transfer for every derivative/window, hence all locally "
         "uniform t->1- limits; prove that the degree-zero endpoint value used by close_left_open_past is the actual left "
         "limit and that all normal jets match the Borel right extension; prove that the actual past residual and every "
         "genuine endpoint normal coefficient vanish outside the official support cylinder so the pointwise support "
         "implication applies globally; derive genuine analytic compact-template derivative bounds for those true jets and "
         "instantiate the landed scale schedule and exact 2^-N derivative-tail budget so the resulting global force is "
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
            "The endpoint Borel/full-spacetime-jet/scale-schedule/traced-residual adapters still depend on caller-supplied endpoint jets or analytic template bounds; close_left_open_past only supplies the prescribed degree-zero endpoint candidate and does not prove the incoming residual converges to it. The endpoint-limit majorant additionally depends on independently proved actual residual derivative bounds with integrable exponent alpha<1. The official-plateau transfer only verifies that such a bound is unchanged on t>=3/4 for T=1 where chi=1 and chi'=0; it does not derive the bound. None of these layers infer endpoint convergence from samples or by themselves prove smooth t=1 force extension; the exact 2^-N omitted-derivative tail budget is conditional on genuine all-order bounds.",
            "The endpoint point-support certificate is an exact implication for a queried point and locally finite active coefficients; it does not prove the unresolved actual residual or all endpoint jets vanish outside the support cylinder, nor does a finite collection of point queries establish global compact support.",
            "The Section 10 compact-support energy adapter is only a fixed-time implication from an independently certified sup bound; it does not prove a uniform bounded-energy estimate as t approaches the blow-up time.",
            "The Section 10 late-origin adapter preserves an upstream-certified origin curl value on 3/4<=t<1 but does not prove that the upstream curl actually blows up.",
            "The Section 7 stress-cone, primary-amplitude ODE and curl-realization modules provide theorem-shaped algebra only; covariance-perturbation and pulse-error-budget inputs remain caller-certified. A generic finite-interval PrimaryODE/Volterra execution adapter with an independent defect diagnostic is present, but its datum/forcing are not yet the actual Proposition 5.5/pulse data and it is not a rigorous numerical certificate. Actual pulse integrals, genuine Lemma-7.4/concentration inputs, localized cylindrical curl hypotheses and a paper-exact wave are still missing.",
            "The Stage-2 cutoff scheduler can preserve theorem-admissible scales beyond binary64 with arbitrary-precision integers and log-space support edges, but it still consumes caller-supplied analytic C[j,m] bounds and constructs only finite prefixes; no infinite schedule, theorem-level local finiteness, or Proposition 5.3 conclusion follows from this representation bridge.",
            "The connected Stage-1 actual-schedule scale chain now uses componentwise field bounds so the current regression resolvent majorant is representable; the original binary64 remainder path still fails closed on overflow, while the wide Decimal adapter carries the same conservative remainderBound/remainderLip algebra to finite values beyond float range and the wide Lambda/C selector carries those values onward without binary64 down-conversion by storing C through its certified exponential exponent. The landed Picard gate additionally certifies the scalar s*remainderBound<=1 and s*remainderLip<=1/2 conditions and exposes the theorem-shaped update/tail bound, but without a genuine coefficient-state backend it still does not materialize the coefficient-space fixed point.",
        ],
        "sources": copy.deepcopy(SOURCE_PINS),
        "stages": copy.deepcopy(STAGES),
    }