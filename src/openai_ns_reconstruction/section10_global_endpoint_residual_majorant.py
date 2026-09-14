"""All-order global small-scale majorant for the actual paper residual.

This is a theorem-shaped Section 10 bridge, not a numerical residual check.  The
pinned ``LocalPaper.Properties.residual_flatness`` statement controls every
natural derivative order and every nonnegative decay power on the inner scaled
region ``radius <= Xext``.  The exterior zero-germ chain recorded by
``section10_exterior_residual_zero_germ`` controls every derivative order on the
complementary paper exterior.  Splitting at the *same* ``Xext`` and taking
``delta = min(delta_inner, qstar)`` therefore yields, for each ``m`` and
``r >= 0``, one global-in-space estimate on the sufficiently small physical-q
sublevel:

    ||D^m R(w)|| <= C * q(w)^r.

The quantifiers over ``m`` and ``r`` remain outside any finite frontier.  This
closes the singular small-scale corner of the late-time residual estimate at
the formal theorem level.  It does not materialize the existential constants,
fields, residual, or force in Python, and it does not by itself cover the
positive-q remainder of the whole ``3/4 <= t < 1`` plateau.
"""
from __future__ import annotations

from dataclasses import dataclass

from openai_ns_reconstruction.section10_exterior_residual_zero_germ import (
    Section10ExteriorResidualZeroGermAdmission,
)
from openai_ns_reconstruction.section10_paper_localization_spine import (
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_REPOSITORY,
    Section10PaperLocalizationSpineAdmission,
)


PINNED_LOCAL_PAPER_FILE = "NavierStokes/LocalPaperTheorem.lean"
PINNED_LOCAL_DOMAIN_FILE = "NavierStokes/LocalPaperDomain.lean"
PINNED_RESIDUAL_FLATNESS_FIELD = "NavierStokes.LocalPaper.Properties.residual_flatness"
PINNED_QSTAR_POS_THEOREM = "NavierStokes.LocalPaperDomain.qstar_pos"
PINNED_OUTER_EDGE = "NavierStokes.LocalPaperDomain.outerEdge"
PINNED_EXTERIOR_CONDITIONS = "NavierStokes.LocalPaperDomain.exterior_conditions"
PINNED_EVIDENCE_KIND = "lean-proof-chain-export"

UNIVERSAL_DERIVATIVE_QUANTIFIER = "forall-m-in-N"
UNIVERSAL_DECAY_QUANTIFIER = "forall-r-in-R-with-0<=r"
GLOBAL_SPATIAL_QUANTIFIER = "forall-w-with-t<1-and-q<delta"
INNER_REGION = "radius(w)<=Xext"
OUTER_REGION = "Xext<=radius(w)"
GLOBAL_DELTA_RULE = "delta=min(delta_inner,qstar)"
MAJORANT_TEMPLATE = "norm(iteratedFDeriv(m,residual,w))<=C*q(w)^r"

REQUIRED_FORMAL_SYMBOLS = (
    PINNED_RESIDUAL_FLATNESS_FIELD,
    PINNED_QSTAR_POS_THEOREM,
    PINNED_OUTER_EDGE,
    PINNED_EXTERIOR_CONDITIONS,
)


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class Section10GlobalEndpointResidualMajorantWitness:
    """Exact metadata for the inner-flatness/exterior-zero global split."""

    application_id: str
    schedule_id: str
    local_velocity_id: str
    local_pressure_id: str
    residual_id: str
    provenance: str

    evidence_kind: str = PINNED_EVIDENCE_KIND
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    local_paper_file: str = PINNED_LOCAL_PAPER_FILE
    local_domain_file: str = PINNED_LOCAL_DOMAIN_FILE
    residual_flatness_field: str = PINNED_RESIDUAL_FLATNESS_FIELD
    qstar_pos_theorem: str = PINNED_QSTAR_POS_THEOREM
    outer_edge_symbol: str = PINNED_OUTER_EDGE
    exterior_conditions_theorem: str = PINNED_EXTERIOR_CONDITIONS
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    derivative_order_quantifier: str = UNIVERSAL_DERIVATIVE_QUANTIFIER
    decay_order_quantifier: str = UNIVERSAL_DECAY_QUANTIFIER
    spatial_quantifier: str = GLOBAL_SPATIAL_QUANTIFIER
    inner_region: str = INNER_REGION
    outer_region: str = OUTER_REGION
    global_delta_rule: str = GLOBAL_DELTA_RULE
    majorant_template: str = MAJORANT_TEMPLATE

    max_derivative_order: int | None = None
    max_decay_order: int | None = None
    sampled_or_fitted_evidence: bool = False
    manufactured_residual: bool = False
    force_defined_as_residual: bool = False

    same_schedule_as_paper_spine_certified: bool = False
    same_local_velocity_as_paper_spine_certified: bool = False
    same_local_pressure_as_paper_spine_certified: bool = False
    same_residual_as_exterior_zero_germ_certified: bool = False
    inner_residual_flatness_at_outer_edge_all_orders_certified: bool = False
    qstar_positive_certified: bool = False
    outer_all_order_zero_jets_certified: bool = False
    inner_outer_radius_split_exhaustive_certified: bool = False
    global_delta_min_positive_certified: bool = False
    global_small_scale_majorant_all_orders_certified: bool = False

    def __post_init__(self) -> None:
        for name in (
            "application_id",
            "schedule_id",
            "local_velocity_id",
            "local_pressure_id",
            "residual_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        expected = {
            "evidence_kind": PINNED_EVIDENCE_KIND,
            "formal_repository": PINNED_FORMAL_REPOSITORY,
            "formal_commit": PINNED_FORMAL_COMMIT,
            "local_paper_file": PINNED_LOCAL_PAPER_FILE,
            "local_domain_file": PINNED_LOCAL_DOMAIN_FILE,
            "residual_flatness_field": PINNED_RESIDUAL_FLATNESS_FIELD,
            "qstar_pos_theorem": PINNED_QSTAR_POS_THEOREM,
            "outer_edge_symbol": PINNED_OUTER_EDGE,
            "exterior_conditions_theorem": PINNED_EXTERIOR_CONDITIONS,
            "dependency_symbols": REQUIRED_FORMAL_SYMBOLS,
            "derivative_order_quantifier": UNIVERSAL_DERIVATIVE_QUANTIFIER,
            "decay_order_quantifier": UNIVERSAL_DECAY_QUANTIFIER,
            "spatial_quantifier": GLOBAL_SPATIAL_QUANTIFIER,
            "inner_region": INNER_REGION,
            "outer_region": OUTER_REGION,
            "global_delta_rule": GLOBAL_DELTA_RULE,
            "majorant_template": MAJORANT_TEMPLATE,
        }
        for name, value in expected.items():
            if getattr(self, name) != value:
                raise ValueError(f"{name} must match the pinned global residual-majorant chain")

        if self.max_derivative_order is not None:
            raise ValueError("finite derivative frontiers cannot certify the all-order majorant")
        if self.max_decay_order is not None:
            raise ValueError("finite decay frontiers cannot certify arbitrary residual decay")
        if self.sampled_or_fitted_evidence is not False:
            raise ValueError("sampled/fitted evidence cannot certify an analytic residual majorant")
        if self.manufactured_residual is not False:
            raise ValueError("manufactured residual evidence is forbidden")
        if self.force_defined_as_residual is not False:
            raise ValueError("f=R is not independent force verification")

        for name in (
            "same_schedule_as_paper_spine_certified",
            "same_local_velocity_as_paper_spine_certified",
            "same_local_pressure_as_paper_spine_certified",
            "same_residual_as_exterior_zero_germ_certified",
            "inner_residual_flatness_at_outer_edge_all_orders_certified",
            "qstar_positive_certified",
            "outer_all_order_zero_jets_certified",
            "inner_outer_radius_split_exhaustive_certified",
            "global_delta_min_positive_certified",
            "global_small_scale_majorant_all_orders_certified",
        ):
            _strict_true(getattr(self, name), name)


@dataclass(frozen=True)
class Section10GlobalEndpointResidualMajorantAdmission:
    """Formal global q-sublevel majorant; all runtime truth flags remain false."""

    spine: Section10PaperLocalizationSpineAdmission
    exterior_zero_germ: Section10ExteriorResidualZeroGermAdmission
    witness: Section10GlobalEndpointResidualMajorantWitness

    def __post_init__(self) -> None:
        if not isinstance(self.spine, Section10PaperLocalizationSpineAdmission):
            raise TypeError("spine must be a Section10PaperLocalizationSpineAdmission")
        if not isinstance(
            self.exterior_zero_germ, Section10ExteriorResidualZeroGermAdmission
        ):
            raise TypeError(
                "exterior_zero_germ must be a Section10ExteriorResidualZeroGermAdmission"
            )
        if not isinstance(self.witness, Section10GlobalEndpointResidualMajorantWitness):
            raise TypeError(
                "witness must be a Section10GlobalEndpointResidualMajorantWitness"
            )
        if self.exterior_zero_germ.spine != self.spine:
            raise ValueError("exterior zero-germ admission must consume the same paper spine")

        parent = self.spine.witness
        exterior = self.exterior_zero_germ.witness
        for child_name, expected in (
            ("schedule_id", parent.schedule_id),
            ("local_velocity_id", parent.local_velocity_id),
            ("local_pressure_id", parent.local_pressure_id),
            ("residual_id", exterior.residual_id),
        ):
            if getattr(self.witness, child_name) != expected:
                raise ValueError(f"{child_name} must be identical to the admitted proof chain")

    @property
    def formal_global_small_scale_residual_majorant_admitted(self) -> bool:
        return True

    @property
    def formal_all_order_residual_decay_family_admitted(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_fields_materialized(self) -> bool:
        return False

    @property
    def actual_majorant_constants_materialized(self) -> bool:
        return False

    @property
    def official_full_late_plateau_majorant_verified(self) -> bool:
        return False

    @property
    def residual_artifact_ready(self) -> bool:
        return False

    @property
    def forcing_artifact_ready(self) -> bool:
        return False

    @property
    def endpoint_residual_closure_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False

    @property
    def full_reconstruction(self) -> bool:
        return False


def admit_global_endpoint_residual_majorant(
    spine: Section10PaperLocalizationSpineAdmission,
    exterior_zero_germ: Section10ExteriorResidualZeroGermAdmission,
    export: Section10GlobalEndpointResidualMajorantWitness,
) -> Section10GlobalEndpointResidualMajorantAdmission:
    """Bind the universal inner/outer split for the same actual paper residual."""

    return Section10GlobalEndpointResidualMajorantAdmission(
        spine=spine,
        exterior_zero_germ=exterior_zero_germ,
        witness=export,
    )
