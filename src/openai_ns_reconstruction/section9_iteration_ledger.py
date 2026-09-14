"""Exact finite-step arithmetic for the pinned Section 9 correction ledger.

This module replays only the stage/exponent arithmetic from the pinned
``NavierStokes/ActualIterationLedger.lean``.  It deliberately does *not*
materialize a wave, mean correction, pressure, residual, or summed local field.

The important boundary is exactness: all stage exponents are represented by
:class:`fractions.Fraction`, and binary floating-point input is rejected.  A
certificate therefore records the paper-selected schedule identities exactly,
rather than inferring cancellation from a small numerical residual.

The finite-step identities certified here are:

* ``sigma(J) = 1/5 + J/10``;
* ``sigma(J+1) - sigma(J) = 1/10``;
* ``gain(h,J) = h*J/10``;
* ``residualWave(J) = J/10 + 7/10``;
* ``residualMean(J) = J/10 + 6/5``;
* both native residual exponents improve by exactly ``1/10`` per cycle; and
* the corresponding physical exponent improves by exactly ``h/10``.

These are finite arithmetic facts.  They are not an existence theorem for the
correction sequence and they do not imply convergence of Eq. (9.21).
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral


PINNED_FORMAL_REPOSITORY = "https://github.com/openai/NavierStokesAndEuler"
PINNED_FORMAL_COMMIT = "f9e8bc5b38b6e212696e8a30e3e91517af887bbd"
PINNED_FORMAL_FILE = "NavierStokes/ActualIterationLedger.lean"
PINNED_SIGMA_SYMBOL = "NavierStokes.ActualIterationLedger.sigma"
PINNED_GAIN_SYMBOL = "NavierStokes.ActualIterationLedger.gain"
PINNED_RESIDUAL_WAVE_SYMBOL = "NavierStokes.ActualIterationLedger.residualWave"
PINNED_RESIDUAL_MEAN_SYMBOL = "NavierStokes.ActualIterationLedger.residualMean"
PINNED_KAPPA = Fraction(1, 100000)


def _stage(value: object, name: str = "stage") -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be a nonnegative integer")
    result = int(value)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


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


def sigma(stage: int) -> Fraction:
    """Pinned accuracy ``1/5 + stage/10`` after ``stage`` cycles."""
    j = _stage(stage)
    return Fraction(1, 5) + Fraction(j, 10)


def input_sigma(positive_stage: int) -> Fraction:
    """Input accuracy for a positive physical increment.

    The formal definition also has a value at zero through natural-number
    subtraction, but that value is explicitly unused by positive-stage
    estimates.  Python rejects zero here so it cannot accidentally be treated
    as a physical increment input.
    """
    j = _stage(positive_stage, "positive_stage")
    if j < 1:
        raise ValueError("positive_stage must be at least 1")
    return sigma(j - 1)


def gain(h: Fraction | int, stage: int) -> Fraction:
    """Pinned common physical gain ``h*stage/10``."""
    hq = _exact_fraction(h, "h")
    j = _stage(stage)
    return hq * Fraction(j, 10)


def residual_wave(stage: int) -> Fraction:
    """Pinned finite-prefix wave residual exponent."""
    j = _stage(stage)
    return Fraction(j, 10) + Fraction(7, 10)


def residual_mean(stage: int) -> Fraction:
    """Pinned finite-prefix mean residual exponent."""
    j = _stage(stage)
    return Fraction(j, 10) + Fraction(6, 5)


def residual_minimum(stage: int) -> Fraction:
    """Exact minimum of the two native residual exponents."""
    wave = residual_wave(stage)
    mean = residual_mean(stage)
    if wave > mean:  # fail closed if the pinned arithmetic is ever changed incorrectly
        raise AssertionError("pinned residual minimum identity no longer holds")
    return wave


@dataclass(frozen=True)
class Section9ResidualImprovementCertificate:
    """Machine-checkable certificate for one finite Section 9 cycle.

    The certificate is intentionally narrow.  It proves only the exact ledger
    improvement from prefix ``J`` to prefix ``J+1``.  It contains no field data
    and therefore cannot be promoted to a residual or convergence certificate.
    """

    stage: int
    h: Fraction
    kappa: Fraction = PINNED_KAPPA

    formal_repository: str = PINNED_FORMAL_REPOSITORY
    formal_commit: str = PINNED_FORMAL_COMMIT
    formal_file: str = PINNED_FORMAL_FILE

    def __post_init__(self) -> None:
        j = _stage(self.stage)
        hq = _exact_fraction(self.h, "h")
        kq = _exact_fraction(self.kappa, "kappa")
        if hq <= 0:
            raise ValueError("h must be strictly positive for residual improvement")
        if kq != PINNED_KAPPA:
            raise ValueError("kappa must equal the pinned paper value 1/100000 exactly")
        if self.formal_repository != PINNED_FORMAL_REPOSITORY:
            raise ValueError("formal_repository must match the pinned source exactly")
        if self.formal_commit != PINNED_FORMAL_COMMIT:
            raise ValueError("formal_commit must match the pinned source exactly")
        if self.formal_file != PINNED_FORMAL_FILE:
            raise ValueError("formal_file must match ActualIterationLedger.lean exactly")
        object.__setattr__(self, "stage", j)
        object.__setattr__(self, "h", hq)
        object.__setattr__(self, "kappa", kq)

        if self.sigma_after - self.sigma_before != Fraction(1, 10):
            raise AssertionError("sigma successor identity failed")
        if self.gain_after - self.gain_before != hq / 10:
            raise AssertionError("physical gain successor identity failed")
        if self.wave_after - self.wave_before != Fraction(1, 10):
            raise AssertionError("wave residual successor identity failed")
        if self.mean_after - self.mean_before != Fraction(1, 10):
            raise AssertionError("mean residual successor identity failed")
        if self.minimum_before != self.wave_before or self.minimum_after != self.wave_after:
            raise AssertionError("residual minimum identity failed")
        if self.physical_minimum_after - self.physical_minimum_before != hq / 10:
            raise AssertionError("physical residual exponent did not improve by exactly h/10")

    @property
    def sigma_before(self) -> Fraction:
        return sigma(self.stage)

    @property
    def sigma_after(self) -> Fraction:
        return sigma(self.stage + 1)

    @property
    def gain_before(self) -> Fraction:
        return gain(self.h, self.stage)

    @property
    def gain_after(self) -> Fraction:
        return gain(self.h, self.stage + 1)

    @property
    def wave_before(self) -> Fraction:
        return residual_wave(self.stage)

    @property
    def wave_after(self) -> Fraction:
        return residual_wave(self.stage + 1)

    @property
    def mean_before(self) -> Fraction:
        return residual_mean(self.stage)

    @property
    def mean_after(self) -> Fraction:
        return residual_mean(self.stage + 1)

    @property
    def minimum_before(self) -> Fraction:
        return residual_minimum(self.stage)

    @property
    def minimum_after(self) -> Fraction:
        return residual_minimum(self.stage + 1)

    @property
    def physical_minimum_before(self) -> Fraction:
        return self.h * self.minimum_before

    @property
    def physical_minimum_after(self) -> Fraction:
        return self.h * self.minimum_after

    @property
    def native_increment(self) -> Fraction:
        return Fraction(1, 10)

    @property
    def physical_increment(self) -> Fraction:
        return self.h / 10

    @property
    def finite_step_only(self) -> bool:
        return True

    @property
    def actual_stage_field_consumed(self) -> bool:
        return False

    @property
    def actual_residual_evaluated(self) -> bool:
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


def certify_residual_improvement(
    stage: int,
    h: Fraction | int,
    *,
    kappa: Fraction | int = PINNED_KAPPA,
) -> Section9ResidualImprovementCertificate:
    """Construct one exact finite-step Section 9 ledger certificate."""
    return Section9ResidualImprovementCertificate(stage=stage, h=h, kappa=kappa)
