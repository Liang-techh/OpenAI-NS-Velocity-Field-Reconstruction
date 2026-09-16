"""Clean-room exact regression for the public profile swirl/pressure identity.

Only short published formula structure is reimplemented here.  It does not
construct the imported leading profile or replay the source checker.
"""
from __future__ import annotations

from fractions import Fraction
from typing import TypeAlias

Exact: TypeAlias = int | Fraction

SOURCE_REPOSITORY = "https://github.com/KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_PUBLIC_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
SOURCE_RECORD_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_PROFILE_ARCHIVE_PATH = "proof_sources/profiles/profiles_body.tex"
SOURCE_PROFILE_SHA256 = "63327f4a6d339d39de096230af1810cb83e7c427ba0b0cc78c576b2750437a95"
SOURCE_CHECKER_PATH = "proof_sources/profiles/exact_checks.py"
SOURCE_CHECKER_SHA256 = "a7ee77024a92d06c1b69055a00b10a7d5fd600ca424fe084e2f4399310cf56e7"
SOURCE_RESULT_PATH = "proof_sources/profiles/exact_checks.json"
SOURCE_RESULT_SHA256 = "5ef2b5f7cee750088143acda8d4a5ee4f0b12f36f0344ab991f5c9333a842e24"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"

PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False
IMPORTED_PROFILE_EXISTENCE_PROVED_HERE = False


def _exact(value: Exact, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


def profile_X_from_cartesian(y1: Exact, y2: Exact) -> Fraction:
    y1 = _exact(y1, name="y1")
    y2 = _exact(y2, name="y2")
    return (y1 * y1 + y2 * y2) / 2


def cartesian_tangential_profile(y1: Exact, y2: Exact, F: Exact) -> tuple[Fraction, Fraction]:
    """Axis-regular tangential factor before the common physical scale."""
    y1 = _exact(y1, name="y1")
    y2 = _exact(y2, name="y2")
    F = _exact(F, name="F")
    return -y2 * F, y1 * F


def squared_norm2(vector: tuple[Exact, Exact]) -> Fraction:
    a = _exact(vector[0], name="vector[0]")
    b = _exact(vector[1], name="vector[1]")
    return a * a + b * b


def regularized_swirl_energy(X: Exact, F: Exact) -> Fraction:
    """Return E^2 after eliminating E=sqrt(2X)F: exactly 2 X F^2."""
    X = _exact(X, name="X")
    F = _exact(F, name="F")
    if X < 0:
        raise ValueError("X must be nonnegative")
    return 2 * X * F * F


def pressure_gradient_from_swirl(F: Exact) -> Fraction:
    """Return the published radial profile derivative Pi_X=F^2."""
    F = _exact(F, name="F")
    return F * F


def swirl_pressure_residual(X: Exact, F: Exact, pi_x: Exact) -> Fraction:
    """Zero iff E^2=2X Pi_X with E eliminated through F."""
    X = _exact(X, name="X")
    pi_x = _exact(pi_x, name="pi_x")
    if X < 0:
        raise ValueError("X must be nonnegative")
    return 2 * X * pi_x - regularized_swirl_energy(X, F)


def cartesian_swirl_energy_residual(y1: Exact, y2: Exact, F: Exact) -> Fraction:
    X = profile_X_from_cartesian(y1, y2)
    vector = cartesian_tangential_profile(y1, y2, F)
    return squared_norm2(vector) - regularized_swirl_energy(X, F)
