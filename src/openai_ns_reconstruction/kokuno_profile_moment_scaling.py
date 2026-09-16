"""Clean-room exact bridge for the public NS-profile five-moment rescaling.

The Kokuno public reader displays the fixed-rectangle substitution ``X = rho*x``
and the induced five-moment scaling

    (M, I, J, S, C_p) =
    (rho*Mhat, rho**(3/2)*Ihat, rho**(3/2)*Jhat, rho*Shat, Chat_p).

This module reimplements only that algebra and a typed binding receipt. It does
not copy the source checker/proof body, construct the imported leading profile,
or infer exterior-field equality from finitely supplied data.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Protocol, TypeAlias, runtime_checkable

Exact: TypeAlias = int | Fraction

SOURCE_REPOSITORY = "https://github.com/KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_PUBLIC_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
SOURCE_PROFILE_ARCHIVE_PATH = "proof_sources/profiles/axis_exact.tex"
SOURCE_PROFILE_SHA256 = "63327f4ab47864ce2508a673aae24873756f491f4597ab8bb177716557937a95"
SOURCE_BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
SOURCE_MANUSCRIPT_SHA256 = "8c8a94add4e7544827292adf2078d1e4d279c226025fdc8c213f9fa9843afa81"

PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False
IMPORTED_PROFILE_EXISTENCE_PROVED_HERE = False


def _exact(value: Exact, *, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction")
    return Fraction(value)


@dataclass(frozen=True)
class FiveMomentState:
    """Exact values of the public five-component matching state at one radius."""

    M: Fraction
    I: Fraction
    J: Fraction
    S: Fraction
    C_p: Fraction

    def __init__(self, *, M: Exact, I: Exact, J: Exact, S: Exact, C_p: Exact) -> None:
        object.__setattr__(self, "M", _exact(M, name="M"))
        object.__setattr__(self, "I", _exact(I, name="I"))
        object.__setattr__(self, "J", _exact(J, name="J"))
        object.__setattr__(self, "S", _exact(S, name="S"))
        object.__setattr__(self, "C_p", _exact(C_p, name="C_p"))

    def residual_against(self, other: "FiveMomentState") -> "FiveMomentState":
        if not isinstance(other, FiveMomentState):
            raise TypeError("other must be a FiveMomentState")
        return FiveMomentState(
            M=self.M - other.M,
            I=self.I - other.I,
            J=self.J - other.J,
            S=self.S - other.S,
            C_p=self.C_p - other.C_p,
        )

    @property
    def is_zero(self) -> bool:
        return all(value == 0 for value in (self.M, self.I, self.J, self.S, self.C_p))


@dataclass(frozen=True)
class FixedRectangleScale:
    """Exact positive scale encoded by ``sqrt_rho`` so ``rho^(3/2)`` stays rational."""

    sqrt_rho: Fraction

    def __init__(self, sqrt_rho: Exact) -> None:
        value = _exact(sqrt_rho, name="sqrt_rho")
        if value <= 0:
            raise ValueError("sqrt_rho must be positive")
        object.__setattr__(self, "sqrt_rho", value)

    @property
    def rho(self) -> Fraction:
        return self.sqrt_rho * self.sqrt_rho

    @property
    def rho_three_halves(self) -> Fraction:
        return self.rho * self.sqrt_rho


def from_fixed_rectangle(hat: FiveMomentState, *, scale: FixedRectangleScale) -> FiveMomentState:
    """Apply the public ``X=rho*x`` five-moment scaling exactly."""

    if not isinstance(hat, FiveMomentState):
        raise TypeError("hat must be a FiveMomentState")
    if not isinstance(scale, FixedRectangleScale):
        raise TypeError("scale must be a FixedRectangleScale")
    return FiveMomentState(
        M=scale.rho * hat.M,
        I=scale.rho_three_halves * hat.I,
        J=scale.rho_three_halves * hat.J,
        S=scale.rho * hat.S,
        C_p=hat.C_p,
    )


def to_fixed_rectangle(state: FiveMomentState, *, scale: FixedRectangleScale) -> FiveMomentState:
    """Invert the public five-moment scaling exactly."""

    if not isinstance(state, FiveMomentState):
        raise TypeError("state must be a FiveMomentState")
    if not isinstance(scale, FixedRectangleScale):
        raise TypeError("scale must be a FixedRectangleScale")
    return FiveMomentState(
        M=state.M / scale.rho,
        I=state.I / scale.rho_three_halves,
        J=state.J / scale.rho_three_halves,
        S=state.S / scale.rho,
        C_p=state.C_p,
    )


@runtime_checkable
class LeadingProfileLike(Protocol):
    """Structural subset implemented by the target ``profiles.LeadingProfile``."""

    name: str
    paper_exact: bool
    provenance: str | None


@dataclass(frozen=True)
class BoundFiveMomentState:
    """Bind exact five-moment data to one target leading-profile object.

    This is an input/provenance receipt only. In particular, a successful bind
    does not prove that the profile exists with the paper's imported properties.
    """

    profile_name: str
    profile_provenance: str
    X_h: Fraction
    moments: FiveMomentState
    data_provenance: str


def bind_profile_moments(
    profile: LeadingProfileLike,
    *,
    X_h: Exact,
    moments: FiveMomentState,
    data_provenance: str,
) -> BoundFiveMomentState:
    """Attach exact matching data to the target ``LeadingProfile`` interface fail-closed."""

    if not isinstance(profile, LeadingProfileLike):
        raise TypeError("profile must implement name, paper_exact, and provenance")
    if not isinstance(profile.name, str) or not profile.name.strip():
        raise ValueError("profile name must be nonempty")
    if not isinstance(profile.paper_exact, bool):
        raise TypeError("profile.paper_exact must be bool")
    if profile.provenance is None or not isinstance(profile.provenance, str) or not profile.provenance.strip():
        raise ValueError("profile must carry explicit provenance")
    radius = _exact(X_h, name="X_h")
    if radius <= 0:
        raise ValueError("matching radius X_h must be positive")
    if not isinstance(moments, FiveMomentState):
        raise TypeError("moments must be a FiveMomentState")
    if not isinstance(data_provenance, str) or not data_provenance.strip():
        raise ValueError("five-moment data require explicit provenance")
    return BoundFiveMomentState(
        profile_name=profile.name.strip(),
        profile_provenance=profile.provenance.strip(),
        X_h=radius,
        moments=moments,
        data_provenance=data_provenance.strip(),
    )
