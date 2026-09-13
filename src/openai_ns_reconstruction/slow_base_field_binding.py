"""Fail-closed slow-box/base-field interface for Sections 6-7.

``slow_labels.freeze_label_phase_data`` already accepts a tangential-base-jet
provider, but that protocol alone does not record *which* slow box and which
upstream base-field object the provider represents. This module closes that
identity seam without inventing the still-missing paper-exact background.

The binding is deliberately sign-free: the ``sigma=+1`` and ``sigma=-1``
labels attached to one Section 6 box must share the same base field. In
addition to the executable binary64 freeze path, the theorem-facing bridge can
bind the same provider identity to the already-landed
``LargeBandPhysicalBaseBinding`` chain ending at ``FinalSlowBase.velocity``.
That second path is metadata/theorem identity only; it never evaluates a deep
large-band field through the binary64 ``DyadicChart`` API.

This is only ``formal-structure``. It does not certify provider values as
paper-exact, replay an upstream field construction, prove Eqs. (7.9)-(7.11),
solve the amplitude ODE, or materialize the paper velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Protocol

from .charts import DyadicChart
from .phase import TangentialBaseJet
from .phase_large_band_physical_base import LargeBandPhysicalBaseBinding
from .slow_labels import (
    ActiveSlowRepresentative,
    FrozenLabelPhaseData,
    freeze_label_phase_data,
)

PRODUCER_KIND = "formal-interface-export"


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _box_key(value: object) -> tuple[int, tuple[int, int, int]]:
    """Normalize a theorem-facing sign-free box key.

    Unlike :class:`SlowLabel`, this metadata identity is allowed to name the
    asymptotic large bands used by the theorem-facing Section 7 adapters.  The
    executable freeze path still receives an ``ActiveSlowRepresentative`` and
    therefore retains the existing binary64 ``ell<=1000`` boundary.
    """

    if not isinstance(value, tuple) or len(value) != 2:
        raise ValueError("box_key must be (ell, a)")
    ell, a = value
    if isinstance(ell, bool) or not isinstance(ell, Integral) or int(ell) < 1:
        raise ValueError("box_key ell must be a positive integer")
    if not isinstance(a, tuple) or len(a) != 3:
        raise ValueError("box_key a must be a length-three integer tuple")
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in a):
        raise ValueError("box_key a must be a length-three integer tuple")
    return int(ell), tuple(int(x) for x in a)


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be interface-export certified true")
    return True


class IdentifiedTangentialBaseJetProvider(Protocol):
    """Tangential-jet provider carrying the opaque identity used by the binding."""

    base_field_provider_id: str

    def tangential_jet(
        self, *, chart: DyadicChart, R: float, Z: float, T: float
    ) -> TangentialBaseJet:
        ...


@dataclass(frozen=True)
class SlowBoxBaseFieldBindingWitness:
    """Opaque identity bridge from one sign-free slow box to an upstream field."""

    box_key: tuple[int, tuple[int, int, int]]
    base_field_family_id: str
    normalized_base_field_id: str
    tangential_jet_provider_id: str
    evaluation_application_id: str
    source_id: str
    source_revision: str
    prepared_instance_id: str
    producer_kind: str
    provenance: str

    box_identity_certified: bool
    sign_independence_certified: bool
    normalized_chart_coordinates_certified: bool
    provider_derivation_certified: bool
    application_certified: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "box_key", _box_key(self.box_key))
        for name in (
            "base_field_family_id",
            "normalized_base_field_id",
            "tangential_jet_provider_id",
            "evaluation_application_id",
            "source_id",
            "source_revision",
            "prepared_instance_id",
            "provenance",
        ):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.producer_kind != PRODUCER_KIND:
            raise ValueError(
                "producer_kind must be formal-interface-export; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )

        for name in (
            "box_identity_certified",
            "sign_independence_certified",
            "normalized_chart_coordinates_certified",
            "provider_derivation_certified",
            "application_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.source_id, self.source_revision, self.prepared_instance_id

    @property
    def binary64_runtime_box_supported(self) -> bool:
        """Whether this theorem identity can also inhabit current ``SlowLabel``."""

        return self.box_key[0] <= 1000


@dataclass(frozen=True)
class TheoremBackedSlowBoxBaseFieldBinding:
    """Bind the provider registry to the landed pinned physical-base theorem path.

    ``LargeBandPhysicalBaseBinding`` already proves, at theorem-metadata level,
    that one admitted slow box uses the ``frequency``/``axial`` fields of the
    selected ``Prepared`` object and that those are the scaled components of the
    same ``FinalSlowBase.velocity``. This composition requires the executable
    provider registry to name that exact box/source/Prepared identity.

    No value is evaluated here. In particular, a deep asymptotic band can pass
    this identity bridge while remaining intentionally unavailable to the
    binary64 ``DyadicChart``/``SlowLabel`` runtime.
    """

    interface: SlowBoxBaseFieldBindingWitness
    physical_base: LargeBandPhysicalBaseBinding

    def __post_init__(self) -> None:
        if not isinstance(self.interface, SlowBoxBaseFieldBindingWitness):
            raise TypeError("interface must be a SlowBoxBaseFieldBindingWitness")
        if not isinstance(self.physical_base, LargeBandPhysicalBaseBinding):
            raise TypeError("physical_base must be a LargeBandPhysicalBaseBinding")
        if self.interface.box_key != self.physical_base.box_key:
            raise ValueError(
                "interface box key must match the theorem physical-base slow box"
            )
        if self.interface.source_key != self.physical_base.source_key:
            raise ValueError(
                "interface source identity must match the theorem physical-base source"
            )
        if self.interface.prepared_key != self.physical_base.prepared_key:
            raise ValueError(
                "interface Prepared identity must match the theorem physical-base source"
            )
        if not self.physical_base.physical_base_identity_theorem_certified:
            raise ValueError("physical-base theorem binding must remain fully certified")

    @property
    def box_key(self) -> tuple[int, tuple[int, int, int]]:
        return self.interface.box_key

    @property
    def source_key(self) -> tuple[str, str]:
        return self.interface.source_key

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.interface.prepared_key

    @property
    def physical_velocity_definition(self) -> str:
        return self.physical_base.witness.physical_velocity_definition

    @property
    def theorem_physical_base_identity_linked(self) -> bool:
        return True

    @property
    def binary64_runtime_box_supported(self) -> bool:
        return self.interface.binary64_runtime_box_supported

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_physical_base_values_materialized(self) -> bool:
        return False

    @property
    def base_jet_values_machine_certified(self) -> bool:
        return False

    @property
    def uniform_phase_estimates_proved(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False


@dataclass(frozen=True)
class BoundFrozenLabelPhaseData:
    """Frozen local phase data plus the sign-free upstream field identity."""

    frozen: FrozenLabelPhaseData
    binding: SlowBoxBaseFieldBindingWitness

    def __post_init__(self) -> None:
        if not isinstance(self.frozen, FrozenLabelPhaseData):
            raise TypeError("frozen must be FrozenLabelPhaseData")
        if not isinstance(self.binding, SlowBoxBaseFieldBindingWitness):
            raise TypeError("binding must be SlowBoxBaseFieldBindingWitness")
        if self.frozen.source.label.box_key != self.binding.box_key:
            raise ValueError("frozen label must use the identical bound slow box")

    @property
    def slow_box_bound_to_base_field_interface(self) -> bool:
        return True

    @property
    def label_signs_share_box_binding(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def base_field_construction_machine_replayed(self) -> bool:
        return False

    @property
    def paper_exact_base_fields_materialized(self) -> bool:
        return False

    @property
    def base_jet_values_machine_certified(self) -> bool:
        return False

    @property
    def uniform_phase_estimates_proved(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False


def freeze_bound_label_phase_data(
    source: ActiveSlowRepresentative,
    provider: IdentifiedTangentialBaseJetProvider,
    binding: SlowBoxBaseFieldBindingWitness,
    *,
    rectangle_radius: float,
    u_star: float,
) -> BoundFrozenLabelPhaseData:
    """Freeze Section 7 data only after slow-box/provider identities agree.

    The returned numerical jet is whatever the caller-supplied provider
    evaluates. Passing this gate therefore establishes interface coherence, not
    paper exactness. The executable path intentionally retains the
    ``ActiveSlowRepresentative``/binary64 band restriction; theorem-facing
    large-band identities are composed separately by
    :class:`TheoremBackedSlowBoxBaseFieldBinding`.
    """

    if not isinstance(source, ActiveSlowRepresentative):
        raise TypeError("source must be ActiveSlowRepresentative")
    if not isinstance(binding, SlowBoxBaseFieldBindingWitness):
        raise TypeError("binding must be SlowBoxBaseFieldBindingWitness")
    if source.label.box_key != binding.box_key:
        raise ValueError("source label must use the identical bound slow box")

    provider_id = getattr(provider, "base_field_provider_id", None)
    if provider_id != binding.tangential_jet_provider_id:
        raise ValueError(
            "provider base_field_provider_id must match tangential_jet_provider_id"
        )

    frozen = freeze_label_phase_data(
        source,
        provider,
        rectangle_radius=rectangle_radius,
        u_star=u_star,
    )
    return BoundFrozenLabelPhaseData(frozen=frozen, binding=binding)
