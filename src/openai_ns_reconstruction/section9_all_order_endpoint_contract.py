"""Fail-closed all-order Section 9 -> Section 10 endpoint contract.

The existing endpoint ladder is intentionally finite.  This module does not
turn a finite prefix into an infinite argument.  Instead it records the exact
*upstream theorem contract* that must be discharged before Section 10 may treat
the residual endpoint data as an all-order input for Taylor--Borel extension.

In particular, admission requires one fixed residual source/revision, exact
endpoint ``T=1``, the official late plateau ``3/4 <= t < 1``, universal
quantification over every natural derivative order, locally-uniform endpoint
limits, an analytic all-order majorant template, endpoint normal-trace matching,
and strict exterior vanishing on one support cylinder.  Finite ``max_order``
claims, sampled/numerical evidence, and manufactured residuals are rejected.

Passing this contract is only ``formal-structure``.  It does not materialize the
paper residual, construct the force, prove the hypotheses, or upgrade any
paper-exact reconstruction flag.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


SECTION10_ENDPOINT_EXACT = Fraction(1, 1)
SECTION10_LATE_PLATEAU_START_EXACT = Fraction(3, 4)
_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem", "lean-formal-export"})
_UNIVERSAL_ORDER_QUANTIFIER = "forall-natural-derivative-orders"


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class Section9AllOrderEndpointResidualWitness:
    """Theorem-facing witness for the missing infinite endpoint step.

    The boolean fields are hypotheses supplied by a theorem/formal export, not
    facts inferred by Python.  ``max_order`` must stay ``None`` and ``finite``
    must stay ``False`` so a finite ladder cannot accidentally satisfy this
    interface.
    """

    source_id: str
    source_revision: str
    residual_identity: str
    evidence_kind: str
    provenance: str
    endpoint: Fraction = SECTION10_ENDPOINT_EXACT
    late_plateau_start: Fraction = SECTION10_LATE_PLATEAU_START_EXACT
    derivative_order_quantifier: str = _UNIVERSAL_ORDER_QUANTIFIER
    analytic_majorant_template: str = ""
    support_cylinder_id: str = ""
    closed_past_all_order_derivative_bounds_certified: bool = False
    locally_uniform_endpoint_limits_all_orders_certified: bool = False
    analytic_majorant_all_orders_certified: bool = False
    endpoint_normal_trace_match_all_orders_certified: bool = False
    support_exterior_zero_neighborhood_all_orders_certified: bool = False
    same_actual_residual_source_certified: bool = False
    closed_past_zero_before_certified: bool = False
    official_late_plateau_certified: bool = False
    infinite_locally_finite_borel_family_certified: bool = False
    max_order: int | None = None
    finite: bool = False

    def __post_init__(self) -> None:
        for name in (
            "source_id",
            "source_revision",
            "residual_identity",
            "provenance",
            "analytic_majorant_template",
            "support_cylinder_id",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem, formal-theorem, or "
                "lean-formal-export; sampled/fitted/numerical evidence is rejected"
            )
        if type(self.endpoint) is not Fraction or self.endpoint != SECTION10_ENDPOINT_EXACT:
            raise ValueError("endpoint must be the exact rational T=1")
        if (
            type(self.late_plateau_start) is not Fraction
            or self.late_plateau_start != SECTION10_LATE_PLATEAU_START_EXACT
        ):
            raise ValueError("late_plateau_start must be the exact rational 3/4")
        if self.derivative_order_quantifier != _UNIVERSAL_ORDER_QUANTIFIER:
            raise ValueError(
                "derivative_order_quantifier must universally quantify all natural orders"
            )
        if self.max_order is not None:
            raise ValueError("max_order must be None; finite endpoint ladders are not all-order evidence")
        if self.finite is not False:
            raise ValueError("finite must be False for an all-order endpoint contract")

        for name in (
            "closed_past_all_order_derivative_bounds_certified",
            "locally_uniform_endpoint_limits_all_orders_certified",
            "analytic_majorant_all_orders_certified",
            "endpoint_normal_trace_match_all_orders_certified",
            "support_exterior_zero_neighborhood_all_orders_certified",
            "same_actual_residual_source_certified",
            "closed_past_zero_before_certified",
            "official_late_plateau_certified",
            "infinite_locally_finite_borel_family_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision


@dataclass(frozen=True)
class Section9AllOrderEndpointContractCertificate:
    """Admitted universal endpoint contract with truth flags kept fail-closed."""

    witness: Section9AllOrderEndpointResidualWitness
    expected_source_key: tuple[str, str]
    expected_residual_identity: str
    status: str = "formal-structure"
    actual_section9_sequence_verified: bool = False
    section9_all_order_endpoint_limits_verified: bool = False
    residual_artifact_ready: bool = False
    forcing_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False
    full_reconstruction: bool = False

    @property
    def formal_all_order_contract_ready(self) -> bool:
        return (
            self.witness.source_key == self.expected_source_key
            and self.witness.residual_identity == self.expected_residual_identity
            and type(self.witness.endpoint) is Fraction
            and self.witness.endpoint == SECTION10_ENDPOINT_EXACT
            and type(self.witness.late_plateau_start) is Fraction
            and self.witness.late_plateau_start == SECTION10_LATE_PLATEAU_START_EXACT
            and self.witness.derivative_order_quantifier == _UNIVERSAL_ORDER_QUANTIFIER
            and self.witness.max_order is None
            and self.witness.finite is False
            and not self.actual_section9_sequence_verified
            and not self.section9_all_order_endpoint_limits_verified
            and not self.residual_artifact_ready
            and not self.forcing_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
            and not self.full_reconstruction
        )


def admit_section9_all_order_endpoint_contract(
    witness: Section9AllOrderEndpointResidualWitness,
    *,
    expected_source_id: str,
    expected_source_revision: str,
    expected_residual_identity: str,
) -> Section9AllOrderEndpointContractCertificate:
    """Admit only a universal formal-source contract for one exact residual.

    ``expected_*`` values are supplied by the already selected upstream residual
    construction.  This exact identity check prevents a formally valid endpoint
    theorem for some other residual revision from being spliced into Section 10.
    """

    if not isinstance(witness, Section9AllOrderEndpointResidualWitness):
        raise TypeError("witness must be a Section9AllOrderEndpointResidualWitness")
    expected_source_key = (
        _nonempty_text(expected_source_id, "expected_source_id"),
        _nonempty_text(expected_source_revision, "expected_source_revision"),
    )
    expected_identity = _nonempty_text(
        expected_residual_identity, "expected_residual_identity"
    )
    if witness.source_key != expected_source_key:
        raise ValueError("all-order endpoint witness source identity mismatch")
    if witness.residual_identity != expected_identity:
        raise ValueError("all-order endpoint witness residual identity mismatch")

    certificate = Section9AllOrderEndpointContractCertificate(
        witness=witness,
        expected_source_key=expected_source_key,
        expected_residual_identity=expected_identity,
    )
    if not certificate.formal_all_order_contract_ready:
        raise ArithmeticError("all-order endpoint contract invariant failed")
    return certificate
