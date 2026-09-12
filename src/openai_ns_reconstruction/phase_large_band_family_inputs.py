"""Fail-closed Section 7 family-input admission on asymptotic slow labels.

The underflow-safe large-band path now has three landed pieces:

* :mod:`phase_large_band_local_base` admits theorem-certified LocalBase bounds;
* :mod:`phase_large_band_base_source` binds those bounds to one stable,
  sign-independent base-source revision; and
* :mod:`phase_large_band_frame` checks the scalar large-band/frame/damping
  inequalities without forming ``Q=2^-ell``.

What was still missing was a typed place to bind an *active signed family
index* and its slow/slot membership to those same source and scalar data.  The
pinned Lean theorem ``BasePhaseGeometry.FamilyData.phase_estimates`` is uniform
on an active carrier/slot only after those membership and family-identity facts
have been proved.  They must not be inferred from a handful of sampled points.

This module supplies that admission gate.  It records only theorem-provenanced
family facts and checks the scalar global hypotheses used on the route from
``FamilyData.phase_estimates`` to ``coordinate_errors``/``damping_error``:
``0<u<=M``, ``r0>0``, ``1/(2 r0)<=M``, ``4 r0 Tg<=M``, the exact source
revision, and exact agreement with the already-landed phase/frame certificate.

No partition, carrier, slot, base field, normal, or oscillatory wave is
constructed here.  Passing this gate is therefore only ``formal-structure``;
it does not establish Eqs. (7.9)--(7.11) for the manuscript construction and
``paper_exact_velocity_available`` remains false.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .phase_large_band_base_source import LargeBandBaseSourceBinding
from .phase_large_band_frame import LargeBandPhaseFrameCertificate
from .phase_large_band_local_base import AsymptoticSlowLabel


_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be finite")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


@dataclass(frozen=True)
class LargeBandFamilyInputWitness:
    """Theorem-facing inputs for one active signed Section 7 family index.

    ``source_id``/``source_revision`` must name the same sign-free base source
    already bound by :class:`LargeBandBaseSourceBinding`.  ``u``, ``B`` and
    ``viscosity`` are repeated deliberately so that downstream scalar geometry
    cannot silently be evaluated with a different parameter set.

    The qualitative flags are theorem/provenance assertions, not numerical
    observations.  In particular, carrier/slot membership is never inferred
    from sampled coordinates.
    """

    label: AsymptoticSlowLabel
    M: float
    u: float
    B: float
    viscosity: float
    reference_radius: float
    Tg: float
    source_id: str
    source_revision: str
    evidence_kind: str
    provenance: str
    family_index_source_identity_certified: bool
    carrier_membership_certified: bool
    slot_membership_certified: bool
    representative_data_identity_certified: bool
    reference_scale_identity_certified: bool
    viscosity_identity_certified: bool

    def __post_init__(self) -> None:
        if not isinstance(self.label, AsymptoticSlowLabel):
            raise TypeError("label must be an AsymptoticSlowLabel")

        M = _finite(self.M, "M")
        u = _finite(self.u, "u")
        B = _finite(self.B, "B")
        viscosity = _finite(self.viscosity, "viscosity")
        r0 = _finite(self.reference_radius, "reference_radius")
        Tg = _finite(self.Tg, "Tg")
        object.__setattr__(self, "M", M)
        object.__setattr__(self, "u", u)
        object.__setattr__(self, "B", B)
        object.__setattr__(self, "viscosity", viscosity)
        object.__setattr__(self, "reference_radius", r0)
        object.__setattr__(self, "Tg", Tg)

        if M < 1.0:
            raise ValueError("M must satisfy M>=1")
        if not 0.0 < u <= M:
            raise ValueError("family geometry requires 0<u<=M")
        if not 0.0 < B <= M:
            raise ValueError("reference scale requires 0<B<=M")
        if not 0.0 <= viscosity <= 4.0:
            raise ValueError("viscosity must lie in [0,4]")
        if r0 <= 0.0:
            raise ValueError("reference_radius must be positive")
        if Tg < 0.0:
            raise ValueError("Tg must be nonnegative")
        if 1.0 / (2.0 * r0) > M:
            raise ValueError("family length hypothesis 1/(2*r0)<=M failed")
        if 4.0 * r0 * Tg > M:
            raise ValueError("family slot hypothesis 4*r0*Tg<=M failed")

        object.__setattr__(self, "source_id", _nonempty_text(self.source_id, "source_id"))
        object.__setattr__(
            self, "source_revision", _nonempty_text(self.source_revision, "source_revision")
        )
        object.__setattr__(
            self, "provenance", _nonempty_text(self.provenance, "theorem provenance")
        )
        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; sampled/fitted evidence is rejected"
            )
        for name in (
            "family_index_source_identity_certified",
            "carrier_membership_certified",
            "slot_membership_certified",
            "representative_data_identity_certified",
            "reference_scale_identity_certified",
            "viscosity_identity_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        return self.label.box_key

    def scalar_family_checks(self) -> dict[str, bool]:
        """Scalar assumptions shared by the pinned family estimates."""

        return {
            "M_ge_1": self.M >= 1.0,
            "u_positive": self.u > 0.0,
            "u_le_M": self.u <= self.M,
            "B_positive": self.B > 0.0,
            "B_le_M": self.B <= self.M,
            "viscosity_in_0_4": 0.0 <= self.viscosity <= 4.0,
            "reference_radius_positive": self.reference_radius > 0.0,
            "length_global_bound": 1.0 / (2.0 * self.reference_radius) <= self.M,
            "slot_global_bound": 4.0 * self.reference_radius * self.Tg <= self.M,
        }


@dataclass(frozen=True)
class LargeBandFamilyInputAdmission:
    """Bind one active signed family witness to source + scalar geometry.

    A successful admission means that one theorem-provenanced signed index has
    the identities/memberships needed to *invoke* the pinned uniform family
    estimates once the cited upstream family theorem exists.  It does not make
    the theorem facts machine-proved and it never promotes one admitted index
    to an all-active-box statement.
    """

    binding: LargeBandBaseSourceBinding
    frame: LargeBandPhaseFrameCertificate
    witness: LargeBandFamilyInputWitness

    def __post_init__(self) -> None:
        if not isinstance(self.binding, LargeBandBaseSourceBinding):
            raise TypeError("binding must be a LargeBandBaseSourceBinding")
        if not isinstance(self.frame, LargeBandPhaseFrameCertificate):
            raise TypeError("frame must be a LargeBandPhaseFrameCertificate")
        if not isinstance(self.witness, LargeBandFamilyInputWitness):
            raise TypeError("witness must be a LargeBandFamilyInputWitness")

        failed = [name for name, ok in self.admission_checks().items() if not ok]
        if failed:
            raise ValueError("uncertified large-band family input admission: " + ", ".join(failed))

    @property
    def label(self) -> AsymptoticSlowLabel:
        return self.witness.label

    @property
    def source_key(self) -> tuple[str, str]:
        return self.witness.source_key

    def admission_checks(self) -> dict[str, bool]:
        w = self.witness
        expected = self.binding.signed_label(w.label.sigma)
        return {
            "signed_label_identity": w.label == expected,
            "box_identity": w.box_key == self.binding.box_key,
            "source_revision_identity": w.source_key == self.binding.source_key,
            "M_identity": w.M == self.binding.source.M == self.frame.scale.M,
            "band_identity": w.label.ell == self.frame.ell == self.binding.admission.ell,
            "u_identity": w.u == self.frame.u,
            "B_identity": w.B == self.frame.B,
            "viscosity_identity": w.viscosity == self.frame.viscosity,
            "base_source_binding_intact": all(self.binding.binding_checks().values()),
            "large_band_scalar_geometry_intact": all(self.frame.scalar_checks().values()),
            "scalar_family_hypotheses": all(w.scalar_family_checks().values()),
            "family_index_source_identity_certified": w.family_index_source_identity_certified is True,
            "carrier_membership_certified": w.carrier_membership_certified is True,
            "slot_membership_certified": w.slot_membership_certified is True,
            "representative_data_identity_certified": w.representative_data_identity_certified is True,
            "reference_scale_identity_certified": w.reference_scale_identity_certified is True,
            "viscosity_identity_certified": w.viscosity_identity_certified is True,
            "theorem_evidence_not_sampled": w.evidence_kind in _ACCEPTED_EVIDENCE,
            "provenance_present": bool(w.provenance),
        }

    @property
    def family_phase_bound(self) -> float:
        """Pinned family-scale normal/slot bound ``phaseConstant(M)/S_*``.

        This is exposed only after the admission checks have passed.  The value
        is still a theorem envelope, not a measured field error.
        """

        return self.frame.family_normal_delta

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_partition_verified(self) -> bool:
        return False

    @property
    def actual_base_fields_verified(self) -> bool:
        return False

    @property
    def family_phase_estimates_machine_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
