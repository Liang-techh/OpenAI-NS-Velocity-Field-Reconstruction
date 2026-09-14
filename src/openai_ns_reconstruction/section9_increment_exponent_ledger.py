"""Exact Section 9 increment exponents and common physical gain.

This module replays the increment-side arithmetic from the pinned
``NavierStokes/ActualIterationLedger.lean``. It is deliberately narrower than
an actual StageEstimates certificate: no wave, mean correction, pressure,
support, residual, or infinite sum is materialized here.

For every positive increment index ``j`` the pinned native exponents are affine
in ``j`` with common successor gain ``1/10``. After multiplying by the paper
parameter ``h`` and applying the fixed component offsets, every component has a
physical exponent at least ``gain(h,j)=h*j/10``. The common gain itself
increases by exactly ``h/10`` from one stage to the next.

All arithmetic is exact ``Fraction`` arithmetic. This is the scale ledger
needed by a later genuine stage-bound/summability consumer; it is not by itself
a proof that the correction sequence exists or that Eq. (9.21) converges.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_LEDGER_FILE = "NavierStokes/ActualIterationLedger.lean"
PINNED_COORDINATE_FILE = "NavierStokes/CoordinateAlgebra.lean"
PINNED_KAPPA = Fraction(1, 100000)

PINNED_SYMBOLS = (
    "NavierStokes.ActualIterationLedger.waveNative_formula",
    "NavierStokes.ActualIterationLedger.wavePressureNative_formula",
    "NavierStokes.ActualIterationLedger.meanNative_formula",
    "NavierStokes.ActualIterationLedger.radialNative_formula",
    "NavierStokes.ActualIterationLedger.fixed_offset_inequalities",
    "NavierStokes.ActualIterationLedger.gain_succ",
    "NavierStokes.CoordinateAlgebra.A",
)

_COMPONENTS = (
    "wavePotential",
    "meanStream",
    "directAngular",
    "wavePressure",
    "meanPressure",
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


def _positive_stage(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError("stage must be a positive integer")
    stage = int(value)
    if stage < 1:
        raise ValueError("stage must be at least 1")
    return stage


def coordinate_A(h: Fraction | int) -> Fraction:
    hq = _exact_fraction(h, "h")
    return Fraction(1, 2) + hq


def gain(h: Fraction | int, stage: int) -> Fraction:
    hq = _exact_fraction(h, "h")
    j = _positive_stage(stage)
    return hq * Fraction(j, 10)


def wave_native(kappa: Fraction | int, stage: int) -> Fraction:
    kq = _exact_fraction(kappa, "kappa")
    j = _positive_stage(stage)
    return Fraction(j, 10) + Fraction(3, 5) - kq


def wave_pressure_native(kappa: Fraction | int, stage: int) -> Fraction:
    kq = _exact_fraction(kappa, "kappa")
    j = _positive_stage(stage)
    return Fraction(j, 10) + Fraction(11, 10) - kq


def mean_native(kappa: Fraction | int, stage: int) -> Fraction:
    kq = _exact_fraction(kappa, "kappa")
    j = _positive_stage(stage)
    return Fraction(j, 10) + Fraction(11, 10) - 2 * kq


def radial_native(kappa: Fraction | int, stage: int) -> Fraction:
    return mean_native(kappa, stage) + 1


@dataclass(frozen=True)
class Section9IncrementExponentCertificate:
    """Exact increment exponent ledger for one positive Section 9 stage."""

    stage: int
    h: Fraction | int
    kappa: Fraction | int = PINNED_KAPPA
    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    dependency_files: tuple[str, str] = (PINNED_LEDGER_FILE, PINNED_COORDINATE_FILE)
    dependency_symbols: tuple[str, ...] = PINNED_SYMBOLS

    def __post_init__(self) -> None:
        stage = _positive_stage(self.stage)
        h = _exact_fraction(self.h, "h")
        kappa = _exact_fraction(self.kappa, "kappa")
        if not Fraction(0) < h < Fraction(1, 2):
            raise ValueError("h must satisfy the paper regime 0 < h < 1/2")
        if kappa != PINNED_KAPPA:
            raise ValueError("kappa must equal the pinned value 1/100000 exactly")
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository does not match the pinned source")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit does not match the pinned source")
        if self.dependency_files != (PINNED_LEDGER_FILE, PINNED_COORDINATE_FILE):
            raise ValueError("dependency_files do not match the pinned increment ledger")
        if self.dependency_symbols != PINNED_SYMBOLS:
            raise ValueError("dependency_symbols do not match the pinned increment ledger")
        object.__setattr__(self, "stage", stage)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "kappa", kappa)

        if self.native_successor_increments != {
            name: Fraction(1, 10) for name in self.native_exponents
        }:
            raise AssertionError("not all pinned native increment exponents advance by 1/10")
        if self.common_gain_successor_increment != h / 10:
            raise AssertionError("common physical gain does not advance by h/10")
        if not self.common_gain_lower_bound_holds:
            raise AssertionError("pinned fixed-offset common-gain inequality failed")
        if min(self.component_slacks.values()) <= 0:
            raise AssertionError("paper-regime component slack must be strictly positive")

    @property
    def native_exponents(self) -> dict[str, Fraction]:
        return {
            "wavePotential": wave_native(self.kappa, self.stage),
            "meanStream": mean_native(self.kappa, self.stage),
            "directAngular": mean_native(self.kappa, self.stage),
            "wavePressure": wave_pressure_native(self.kappa, self.stage),
            "meanPressure": mean_native(self.kappa, self.stage),
            "radialMean": radial_native(self.kappa, self.stage),
        }

    @property
    def successor_native_exponents(self) -> dict[str, Fraction]:
        return {
            "wavePotential": wave_native(self.kappa, self.stage + 1),
            "meanStream": mean_native(self.kappa, self.stage + 1),
            "directAngular": mean_native(self.kappa, self.stage + 1),
            "wavePressure": wave_pressure_native(self.kappa, self.stage + 1),
            "meanPressure": mean_native(self.kappa, self.stage + 1),
            "radialMean": radial_native(self.kappa, self.stage + 1),
        }

    @property
    def native_successor_increments(self) -> dict[str, Fraction]:
        after = self.successor_native_exponents
        return {name: after[name] - value for name, value in self.native_exponents.items()}

    @property
    def fixed_offsets(self) -> dict[str, Fraction]:
        return {
            "wavePotential": self.h,
            "meanStream": Fraction(0),
            "directAngular": Fraction(0),
            "wavePressure": 2 * coordinate_A(self.h),
            "meanPressure": Fraction(0),
        }

    @property
    def physical_component_exponents(self) -> dict[str, Fraction]:
        """Replay the five expressions on the RHS of fixed_offset_inequalities."""
        native = self.native_exponents
        offsets = self.fixed_offsets
        return {
            "wavePotential": self.h * native["wavePotential"] - self.h + offsets["wavePotential"],
            "meanStream": self.h * native["meanStream"] + offsets["meanStream"],
            "directAngular": self.h * native["directAngular"] + offsets["directAngular"],
            "wavePressure": (
                self.h * native["wavePressure"]
                - 2 * coordinate_A(self.h)
                + offsets["wavePressure"]
            ),
            "meanPressure": self.h * native["meanPressure"] + offsets["meanPressure"],
        }

    @property
    def common_gain(self) -> Fraction:
        return gain(self.h, self.stage)

    @property
    def successor_common_gain(self) -> Fraction:
        return gain(self.h, self.stage + 1)

    @property
    def common_gain_successor_increment(self) -> Fraction:
        return self.successor_common_gain - self.common_gain

    @property
    def component_slacks(self) -> dict[str, Fraction]:
        return {
            name: exponent - self.common_gain
            for name, exponent in self.physical_component_exponents.items()
        }

    @property
    def common_gain_lower_bound_holds(self) -> bool:
        return set(self.physical_component_exponents) == set(_COMPONENTS) and all(
            self.common_gain <= exponent
            for exponent in self.physical_component_exponents.values()
        )

    @property
    def exact_increment_ledger_machine_checked(self) -> bool:
        return True

    @property
    def finite_stage_only(self) -> bool:
        return True

    @property
    def actual_stage_field_consumed(self) -> bool:
        return False

    @property
    def support_estimates_consumed(self) -> bool:
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
