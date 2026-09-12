"""Fail-closed base-source identity binding for large-band Section 6/7 labels.

``slow_labels.TangentialBaseJetProvider`` is deliberately sign-free: the base
jet is queried from ``(chart,R,Z,T)`` and the wave sign ``sigma`` is introduced
only after that jet has been frozen.  The underflow-safe
:mod:`phase_large_band_local_base` layer, however, previously stopped at an
``AsymptoticSlowLabel`` and a theorem-level ``LocalBase`` admission.  A future
paper-exact provider therefore had no typed place to state that the two labels
``(ell,a,-1)`` and ``(ell,a,+1)`` use the *same* base-field source.

This module supplies that missing identity gate.  It binds one already-admitted
large-band LocalBase witness to a stable source/revision key and exposes the two
sign labels from the common sign-free box key.  The gate rejects sampled or
fitted evidence and rejects any attempt to change the band, box, ``M`` or base
source across the sign pair.

No field values are evaluated here, and the source assertions are not proved by
this Python module.  Until a genuine Proposition 5.5/background theorem witness
is connected, this remains strictly ``formal-structure`` and
``paper_exact_velocity_available`` is false.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral

from .phase_large_band_local_base import (
    AsymptoticSlowLabel,
    LargeBandLocalBaseAdmission,
)


_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _index3(value: object) -> tuple[int, int, int]:
    if not isinstance(value, tuple) or len(value) != 3:
        raise ValueError("a must be a length-three integer tuple")
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in value):
        raise ValueError("a must be a length-three integer tuple")
    return tuple(int(x) for x in value)


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _finite_M(value: object) -> float:
    if isinstance(value, bool):
        raise ValueError("M must be finite with M>=1")
    out = float(value)
    if not math.isfinite(out) or out < 1.0:
        raise ValueError("M must be finite with M>=1")
    return out


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class LargeBandBaseSourceWitness:
    """Theorem-facing identity of the base source for one sign-free slow box.

    ``source_id`` names the provider/family and ``source_revision`` freezes the
    exact theorem/data revision.  The three qualitative facts are separate on
    purpose: a caller must certify that the source obeys the normalized base
    provider contract, is independent of the wave sign, and is exactly the
    source to which the admitted LocalBase theorem applies.

    This is provenance metadata, not a machine proof of the cited theorem.
    """

    ell: int
    a: tuple[int, int, int]
    M: float
    source_id: str
    source_revision: str
    evidence_kind: str
    provenance: str
    provider_contract_certified: bool
    sign_independent_source_certified: bool
    local_base_source_identity_certified: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "ell", _positive_integer(self.ell, "ell"))
        object.__setattr__(self, "a", _index3(self.a))
        object.__setattr__(self, "M", _finite_M(self.M))
        object.__setattr__(self, "source_id", _nonempty_text(self.source_id, "source_id"))
        object.__setattr__(
            self, "source_revision", _nonempty_text(self.source_revision, "source_revision")
        )
        object.__setattr__(self, "provenance", _nonempty_text(self.provenance, "theorem provenance"))
        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; sampled/fitted evidence is rejected"
            )
        for name in (
            "provider_contract_certified",
            "sign_independent_source_certified",
            "local_base_source_identity_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        return self.ell, self.a

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision


@dataclass(frozen=True)
class LargeBandBaseSourceBinding:
    """Bind one LocalBase admission to one sign-independent base source.

    The LocalBase witness may itself carry either sign because Section 6 labels
    are signed.  This binding intentionally discards that sign when identifying
    the base source, then recreates the two allowed labels from the common
    ``box_key``.  Hence downstream wave code cannot silently attach the two
    signs of one slow box to unrelated base-field revisions.
    """

    admission: LargeBandLocalBaseAdmission
    source: LargeBandBaseSourceWitness

    def __post_init__(self) -> None:
        if not isinstance(self.admission, LargeBandLocalBaseAdmission):
            raise TypeError("admission must be a LargeBandLocalBaseAdmission")
        if not isinstance(self.source, LargeBandBaseSourceWitness):
            raise TypeError("source must be a LargeBandBaseSourceWitness")
        label = self.admission.witness.label
        if self.source.box_key != label.box_key:
            raise ValueError("base source box key must match the admitted LocalBase slow box")
        if self.source.M != self.admission.witness.M:
            raise ValueError("base source M must match the admitted LocalBase M exactly")
        failed = [name for name, ok in self.binding_checks().items() if not ok]
        if failed:
            raise ValueError("uncertified large-band base-source binding: " + ", ".join(failed))

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        return self.source.box_key

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source.source_key

    @property
    def sign_pair(self) -> tuple[AsymptoticSlowLabel, AsymptoticSlowLabel]:
        ell, a = self.box_key
        return (
            AsymptoticSlowLabel(ell=ell, a=a, sigma=-1),
            AsymptoticSlowLabel(ell=ell, a=a, sigma=1),
        )

    def signed_label(self, sigma: int) -> AsymptoticSlowLabel:
        if sigma not in (-1, 1):
            raise ValueError("sigma must be +1 or -1")
        ell, a = self.box_key
        return AsymptoticSlowLabel(ell=ell, a=a, sigma=sigma)

    def binding_checks(self) -> dict[str, bool]:
        label = self.admission.witness.label
        return {
            "box_identity": self.source.box_key == label.box_key,
            "band_identity": self.source.ell == self.admission.ell,
            "M_identity": self.source.M == self.admission.witness.M,
            "localbase_admission_intact": all(self.admission.admission_checks().values()),
            "provider_contract_certified": self.source.provider_contract_certified is True,
            "sign_independent_source_certified": self.source.sign_independent_source_certified is True,
            "local_base_source_identity_certified": self.source.local_base_source_identity_certified is True,
            "theorem_evidence_not_sampled": self.source.evidence_kind in _ACCEPTED_EVIDENCE,
            "stable_source_identity_present": bool(self.source.source_id and self.source.source_revision),
            "provenance_present": bool(self.source.provenance),
        }

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_base_fields_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
