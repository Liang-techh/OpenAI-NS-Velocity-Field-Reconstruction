"""Executable finite-prefix field assembly for the Section 9 iteration.

This module materializes one piece of the finite-stage ledger that was
previously represented only by opaque formal identities.  The pinned formal
construction defines

``finiteVelocity J = curl (sum_{j=0}^J potentialStages j)
                    + sum_{j=0}^J directStages j``

and

``finitePressure J = sum_{j=0}^J pressureStages j``.

The equalities are the literal ``MixedDiagonalResidual.uncutVelocity`` and
``DiagonalJetBounds.uncutPrefix`` definitions used by
``WholeDomainStageBounds.finiteVelocity/finitePressure``.  Spatial curl is
linear, so a materialized stage may expose its already-computed curl of the
potential together with the direct velocity and pressure fields.

This is an execution layer, not a source of missing paper data.  It accepts
only a contiguous stage family whose stable family identifiers agree with an
already-admitted :class:`Section9ActualStageEstimatesAdmission`.  It never
fills a missing stage with zero and never infers a field from samples.  The
current repository still lacks the actual materialized stage callables, so
this constructor does not make the reconstruction paper-exact by itself.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral
from typing import Callable, Generic, Iterable, TypeVar

from .section9_actual_stage_estimates import Section9ActualStageEstimatesAdmission


PointT = TypeVar("PointT")
Vector3 = tuple[float, float, float]
VectorField = Callable[[PointT], Vector3]
ScalarField = Callable[[PointT], float]

PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_WHOLE_DOMAIN_FILE = "NavierStokes/WholeDomainActualStageBounds.lean"
PINNED_MIXED_FILE = "NavierStokes/MixedDiagonalResidual.lean"
PINNED_PREFIX_FILE = "NavierStokes/DiagonalJetBounds.lean"
PINNED_FINITE_VELOCITY = "NavierStokes.WholeDomainStageBounds.finiteVelocity"
PINNED_FINITE_PRESSURE = "NavierStokes.WholeDomainStageBounds.finitePressure"
PINNED_UNCUT_VELOCITY = "NavierStokes.MixedDiagonalResidual.uncutVelocity"
PINNED_UNCUT_PREFIX = "NavierStokes.DiagonalJetBounds.uncutPrefix"
MATERIALIZATION_KIND = "machine-materialized-stage-field"


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _finite_scalar(value: object, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must evaluate to a finite real")
    return value


def _finite_vector(value: object, name: str) -> Vector3:
    if not isinstance(value, (tuple, list)) or len(value) != 3:
        raise ValueError(f"{name} must evaluate to exactly three components")
    return tuple(_finite_scalar(v, f"{name}[{i}]") for i, v in enumerate(value))  # type: ignore[return-value]


def _add3(a: Vector3, b: Vector3) -> Vector3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


@dataclass(frozen=True)
class Section9MaterializedStage(Generic[PointT]):
    """One genuinely supplied stage field, with no reconstruction by this module.

    ``curl_potential`` is the spatial curl of the stage potential, not the
    potential itself.  This keeps the runtime representation faithful to the
    exact finite-velocity formula without inventing a numerical curl operator.
    The actual stage producer remains responsible for constructing that field.
    """

    index: int
    potential_family_id: str
    direct_family_id: str
    pressure_family_id: str
    stage_id: str
    curl_potential: VectorField[PointT]
    direct_velocity: VectorField[PointT]
    pressure: ScalarField[PointT]
    materialization_kind: str
    provenance: str
    field_materialization_certified: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "index", _natural(self.index, "index"))
        for name in (
            "potential_family_id",
            "direct_family_id",
            "pressure_family_id",
            "stage_id",
            "provenance",
        ):
            object.__setattr__(self, name, _text(getattr(self, name), name))
        if self.materialization_kind != MATERIALIZATION_KIND:
            raise ValueError(
                "materialization_kind must identify a machine-materialized stage field; "
                "sampled/fitted/toy field payloads are rejected"
            )
        if self.field_materialization_certified is not True:
            raise ValueError("field_materialization_certified must be exactly true")
        if not callable(self.curl_potential) or not callable(self.direct_velocity) or not callable(self.pressure):
            raise TypeError("stage field payloads must be callable fields")

    def curl_value(self, point: PointT) -> Vector3:
        return _finite_vector(self.curl_potential(point), "curl_potential")

    def direct_value(self, point: PointT) -> Vector3:
        return _finite_vector(self.direct_velocity(point), "direct_velocity")

    def pressure_value(self, point: PointT) -> float:
        return _finite_scalar(self.pressure(point), "pressure")


@dataclass(frozen=True)
class Section9FinitePrefixFields(Generic[PointT]):
    """Literal stages ``0,...,J`` assembled by the pinned finite-prefix algebra."""

    estimates: Section9ActualStageEstimatesAdmission
    stages: tuple[Section9MaterializedStage[PointT], ...]
    J: int

    def __post_init__(self) -> None:
        if not isinstance(self.estimates, Section9ActualStageEstimatesAdmission):
            raise TypeError("estimates must be a Section9ActualStageEstimatesAdmission")
        object.__setattr__(self, "J", _natural(self.J, "J"))
        stages = tuple(self.stages)
        object.__setattr__(self, "stages", stages)
        if len(stages) != self.J + 1:
            raise ValueError("finite prefix must contain exactly the stages 0,...,J")
        expected = tuple(range(self.J + 1))
        actual = tuple(stage.index for stage in stages)
        if actual != expected:
            raise ValueError("stage indices must be contiguous and ordered exactly as 0,...,J")
        if len({stage.stage_id for stage in stages}) != len(stages):
            raise ValueError("every materialized stage must have a distinct stage_id")

        witness = self.estimates.witness
        for stage in stages:
            if stage.potential_family_id != witness.finite_velocity_family_id:
                raise ValueError("potential stage family is not the admitted finite-velocity family")
            if stage.direct_family_id != witness.finite_velocity_family_id:
                raise ValueError("direct stage family is not the admitted finite-velocity family")
            if stage.pressure_family_id != witness.finite_pressure_family_id:
                raise ValueError("pressure stage family is not the admitted finite-pressure family")

    @classmethod
    def from_stages(
        cls,
        estimates: Section9ActualStageEstimatesAdmission,
        stages: Iterable[Section9MaterializedStage[PointT]],
    ) -> "Section9FinitePrefixFields[PointT]":
        materialized = tuple(stages)
        if not materialized:
            raise ValueError("a finite prefix must include stage zero")
        return cls(estimates=estimates, stages=materialized, J=len(materialized) - 1)

    def velocity(self, point: PointT) -> Vector3:
        """Evaluate ``curl(sum potential_j) + sum direct_j`` for ``j=0..J``."""
        curl_sum: Vector3 = (0.0, 0.0, 0.0)
        direct_sum: Vector3 = (0.0, 0.0, 0.0)
        for stage in self.stages:
            curl_sum = _add3(curl_sum, stage.curl_value(point))
            direct_sum = _add3(direct_sum, stage.direct_value(point))
        return _add3(curl_sum, direct_sum)

    def pressure(self, point: PointT) -> float:
        """Evaluate the literal pressure prefix ``sum_{j=0}^J pressure_j``."""
        value = 0.0
        for stage in self.stages:
            value += stage.pressure_value(point)
        if not math.isfinite(value):
            raise ValueError("finite pressure prefix overflowed to a nonfinite value")
        return value

    @property
    def stage_count(self) -> int:
        return self.J + 1

    @property
    def exact_prefix_shape_materialized(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_paper_stage_fields_available(self) -> bool:
        # This class checks supplied materialized fields but does not manufacture
        # the missing paper stages or certify their mathematical provenance.
        return False

    @property
    def section9_infinite_iteration_closed(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
