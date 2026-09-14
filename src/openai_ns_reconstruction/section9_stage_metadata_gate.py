"""Fail-closed finite-stage replay of Section 9 ``StageMetadata`` premises.

The pinned formal ``ActualIterationLedger.StageMetadata`` structure does not
construct a correction stage.  Instead, for every positive stage it requires
five native exponent lower bounds plus the two wave/pressure shift lower
bounds.  Those hypotheses are exactly what allows
``StageMetadata.gain_inequalities`` to turn the paper-selected native exponents
into the common physical gain used by the Section 9 smooth-sum estimates.

This module closes only that arithmetic/identity seam.  It consumes an already
materialized :class:`Section9MaterializedStage` together with source-revision-
bound exact metadata for that same stage, checks the seven formal inequalities
with ``Fraction`` arithmetic, and independently rechecks the five resulting
common-gain inequalities.  It never derives exponent metadata from sampled
field values, never fills a missing stage with zero, and never claims that the
repository currently contains the paper's actual correction stages.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Generic, TypeVar

from .section9_finite_prefix_fields import Section9MaterializedStage
from .section9_increment_exponent_ledger import (
    PINNED_FORMAL_COMMIT,
    PINNED_FORMAL_REPOSITORY,
    PINNED_KAPPA,
    Section9IncrementExponentCertificate,
    coordinate_A,
    mean_native,
    wave_native,
    wave_pressure_native,
)


PointT = TypeVar("PointT")
METADATA_KIND = "machine-materialized-exact-stage-metadata"
PINNED_FORMAL_FILE = "NavierStokes/ActualIterationLedger.lean"
PINNED_STAGE_METADATA = "NavierStokes.ActualIterationLedger.StageMetadata"
PINNED_GAIN_INEQUALITIES = (
    "NavierStokes.ActualIterationLedger.StageMetadata.gain_inequalities"
)
REQUIRED_FORMAL_SYMBOLS = (
    PINNED_STAGE_METADATA,
    PINNED_GAIN_INEQUALITIES,
    "NavierStokes.ActualIterationLedger.waveNative_formula",
    "NavierStokes.ActualIterationLedger.wavePressureNative_formula",
    "NavierStokes.ActualIterationLedger.meanNative_formula",
    "NavierStokes.ActualIterationLedger.fixed_offset_inequalities",
)


def _exact_fraction(value: object, name: str) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact integer/Fraction, not bool")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Fraction(int(value), 1)
    raise TypeError(
        f"{name} must be an exact integer/Fraction; float/Decimal approximations are rejected"
    )


def _positive_stage(value: object, name: str = "stage_index") -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be a positive integer")
    result = int(value)
    if result < 1:
        raise ValueError(f"{name} must be at least 1 because StageMetadata applies to increments")
    return result


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _git_sha(value: object) -> str:
    revision = _nonempty_text(value, "source_revision")
    if len(revision) != 40 or any(ch not in "0123456789abcdef" for ch in revision):
        raise ValueError("source_revision must be an exact 40-character lowercase git SHA")
    return revision


@dataclass(frozen=True)
class Section9ExactStageMetadata:
    """Exact metadata exported for one already-produced positive stage.

    Stable stage/family identities are repeated intentionally.  The admission
    object cross-checks them against the supplied materialized stage so metadata
    from one stage cannot silently certify another.
    """

    stage_index: int
    stage_id: str
    potential_family_id: str
    direct_family_id: str
    pressure_family_id: str
    source_revision: str
    producer_artifact_id: str
    provenance: str

    wave_potential_alpha: Fraction | int
    wave_potential_shift: Fraction | int
    mean_stream_alpha: Fraction | int
    direct_angular_alpha: Fraction | int
    wave_pressure_alpha: Fraction | int
    wave_pressure_shift: Fraction | int
    mean_pressure_alpha: Fraction | int

    metadata_kind: str = METADATA_KIND
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE
    dependency_symbols: tuple[str, ...] = REQUIRED_FORMAL_SYMBOLS

    def __post_init__(self) -> None:
        object.__setattr__(self, "stage_index", _positive_stage(self.stage_index))
        for name in (
            "stage_id",
            "potential_family_id",
            "direct_family_id",
            "pressure_family_id",
            "producer_artifact_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))
        object.__setattr__(self, "source_revision", _git_sha(self.source_revision))
        for name in (
            "wave_potential_alpha",
            "wave_potential_shift",
            "mean_stream_alpha",
            "direct_angular_alpha",
            "wave_pressure_alpha",
            "wave_pressure_shift",
            "mean_pressure_alpha",
        ):
            object.__setattr__(self, name, _exact_fraction(getattr(self, name), name))

        if self.metadata_kind != METADATA_KIND:
            raise ValueError(
                "metadata_kind must identify machine-materialized exact stage metadata; "
                "sampled/fitted/toy metadata are rejected"
            )
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository does not match the pinned source")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit does not match the pinned source")
        if self.formal_file != PINNED_FORMAL_FILE:
            raise ValueError("formal_file must be ActualIterationLedger.lean")
        if self.dependency_symbols != REQUIRED_FORMAL_SYMBOLS:
            raise ValueError("dependency_symbols do not match the pinned StageMetadata chain")


@dataclass(frozen=True)
class Section9StageMetadataAdmission(Generic[PointT]):
    """Exact admission of the seven finite-stage ``StageMetadata`` premises."""

    stage: Section9MaterializedStage[PointT]
    metadata: Section9ExactStageMetadata
    h: Fraction | int
    kappa: Fraction | int = PINNED_KAPPA

    def __post_init__(self) -> None:
        if not isinstance(self.stage, Section9MaterializedStage):
            raise TypeError("stage must be a Section9MaterializedStage")
        if not isinstance(self.metadata, Section9ExactStageMetadata):
            raise TypeError("metadata must be Section9ExactStageMetadata")

        h = _exact_fraction(self.h, "h")
        kappa = _exact_fraction(self.kappa, "kappa")
        if not Fraction(0) < h < Fraction(1, 2):
            raise ValueError("h must satisfy the paper regime 0 < h < 1/2")
        if kappa != PINNED_KAPPA:
            raise ValueError("kappa must equal the pinned value 1/100000 exactly")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "kappa", kappa)

        md = self.metadata
        if self.stage.index != md.stage_index:
            raise ValueError("stage index is cross-wired against exact stage metadata")
        if self.stage.index < 1:
            raise ValueError("StageMetadata admission requires a positive increment stage")
        if self.stage.stage_id != md.stage_id:
            raise ValueError("stage_id is cross-wired against exact stage metadata")
        if self.stage.potential_family_id != md.potential_family_id:
            raise ValueError("potential family is cross-wired against exact stage metadata")
        if self.stage.direct_family_id != md.direct_family_id:
            raise ValueError("direct family is cross-wired against exact stage metadata")
        if self.stage.pressure_family_id != md.pressure_family_id:
            raise ValueError("pressure family is cross-wired against exact stage metadata")

        violations = [
            ("wave_potential_alpha", md.wave_potential_alpha, wave_native(kappa, md.stage_index)),
            ("mean_stream_alpha", md.mean_stream_alpha, mean_native(kappa, md.stage_index)),
            ("direct_angular_alpha", md.direct_angular_alpha, mean_native(kappa, md.stage_index)),
            ("wave_pressure_alpha", md.wave_pressure_alpha, wave_pressure_native(kappa, md.stage_index)),
            ("mean_pressure_alpha", md.mean_pressure_alpha, mean_native(kappa, md.stage_index)),
            ("wave_potential_shift", md.wave_potential_shift, -h),
            ("wave_pressure_shift", md.wave_pressure_shift, -2 * coordinate_A(h)),
        ]
        for name, actual, required in violations:
            if actual < required:
                raise ValueError(
                    f"{name} violates pinned StageMetadata lower bound: "
                    f"actual={actual}, required>={required}"
                )

        if not all(
            self.ledger.common_gain <= exponent
            for exponent in self.achieved_physical_exponents.values()
        ):
            raise AssertionError("StageMetadata premises did not imply all five pinned gain inequalities")

    @property
    def ledger(self) -> Section9IncrementExponentCertificate:
        return Section9IncrementExponentCertificate(
            stage=self.metadata.stage_index,
            h=self.h,
            kappa=self.kappa,
        )

    @property
    def achieved_physical_exponents(self) -> dict[str, Fraction]:
        md = self.metadata
        A = coordinate_A(self.h)
        return {
            "wavePotential": self.h * md.wave_potential_alpha + md.wave_potential_shift + self.h,
            "meanStream": self.h * md.mean_stream_alpha,
            "directAngular": self.h * md.direct_angular_alpha,
            "wavePressure": self.h * md.wave_pressure_alpha + md.wave_pressure_shift + 2 * A,
            "meanPressure": self.h * md.mean_pressure_alpha,
        }

    @property
    def gain_slacks(self) -> dict[str, Fraction]:
        g = self.ledger.common_gain
        return {name: exponent - g for name, exponent in self.achieved_physical_exponents.items()}

    @property
    def exact_stage_metadata_premises_verified(self) -> bool:
        return True

    @property
    def materialized_stage_identity_consumed(self) -> bool:
        return True

    @property
    def finite_positive_stage_only(self) -> bool:
        return True

    @property
    def field_values_or_norms_evaluated(self) -> bool:
        return False

    @property
    def support_estimates_consumed(self) -> bool:
        return False

    @property
    def actual_paper_stage_certified(self) -> bool:
        return False

    @property
    def infinite_correction_sequence_certified(self) -> bool:
        return False

    @property
    def eq_9_21_summed_field_certified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
