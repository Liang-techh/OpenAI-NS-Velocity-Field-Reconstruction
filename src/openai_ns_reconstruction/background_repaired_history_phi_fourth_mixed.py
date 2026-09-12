"""Hierarchy-owned repaired fourth-mixed ``phi_n`` jets for Section 5.

The landed compact Lemma 5.2 adapter now supplies the three total-order-four
``phi_n`` derivatives required by the analytic second-eta preceding-diffusion
row.  This module makes that stronger angular layer authoritative in the same
coefficient source that already owns repaired fifth-mixed ``U_n`` data.

The fourth-mixed phi provider is authoritative: its third- and second-mixed
projections are the inherited phi providers, while the fifth-mixed U provider
continues to own all lower U projections and all beta data are still derived
only through the landed analytic Eq. (5.2) adapters.  Missing genuine order-zero
strong data therefore fail closed instead of being replaced by a caller table.

This remains Stage-2 ``formal-structure`` infrastructure.  The unrepaired base
jets, moment/patch eta jets, normalization C, and the genuine leading profile
remain upstream inputs and are not promoted to paper-exact data here.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .background_moment_repair_fifth_mixed_jets import (
    BaseFifthMixedJetProvider,
    Lemma52RepairedFifthMixedJetAdapter,
)
from .background_moment_repair_phi_fourth_mixed_jets import (
    BasePhiFourthMixedJetProvider,
    Lemma52RepairedPhiFourthMixedJetAdapter,
)
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_preceding_diffusion_parameter_jet import ProfileThirdMixedJet
from .background_preceding_diffusion_second_parameter_jet import (
    ProfileFourthMixedJet,
)
from .background_regular_flux_fourth_mixed_jets import AxialFifthMixedJet
from .background_regular_flux_second_jets import AxialThirdMixedJet
from .background_regular_flux_third_mixed_jets import AxialFourthMixedJet
from .background_repaired_history_fifth_mixed import (
    Section5FifthMixedCoefficientJetSource,
    Section5LowerHistoryFifthMixedHierarchy,
)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


@dataclass(frozen=True)
class Section5PhiFourthMixedCoefficientJetSource(
    Section5FifthMixedCoefficientJetSource
):
    """Coefficient source with one authoritative repaired fourth-mixed phi jet."""

    phi_fourth_mixed_jet_provider: BasePhiFourthMixedJetProvider | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if not callable(self.phi_fourth_mixed_jet_provider):
            raise TypeError("phi_fourth_mixed_jet_provider must be callable")

    @classmethod
    def from_lemma52_repair_full_strong_mixed(
        cls,
        order: int,
        *,
        C: float,
        profile: Lemma52RepairedProfileAdapter,
        base_phi_fourth_mixed_jet: BasePhiFourthMixedJetProvider,
        base_U_fifth_mixed_jet: BaseFifthMixedJetProvider,
        provenance: str | None = None,
    ) -> "Section5PhiFourthMixedCoefficientJetSource":
        """Own repaired phi-fourth and U-fifth jets for one positive order.

        Lower angular and axial providers are exact projections of those two
        authoritative strong providers.  No independent beta derivative table
        is accepted; inherited hierarchy methods continue to derive beta only
        through analytic Eq. (5.2).
        """

        order = _positive_order(order)
        phi = Lemma52RepairedPhiFourthMixedJetAdapter(
            profile,
            C,
            base_phi_fourth_mixed_jet,
        )
        axial = Lemma52RepairedFifthMixedJetAdapter(
            profile,
            base_U_fifth_mixed_jet,
        )

        def phi_fourth_provider(X: float, eta: float) -> ProfileFourthMixedJet:
            return phi.phi_fourth_mixed_jet(X, eta)

        def phi_third_provider(X: float, eta: float) -> ProfileThirdMixedJet:
            return phi_fourth_provider(X, eta).third()

        def phi_second_provider(X: float, eta: float):
            return phi_third_provider(X, eta).second()

        def fifth_provider(X: float, eta: float) -> AxialFifthMixedJet:
            return axial.U_fifth_mixed_jet(X, eta)

        def fourth_provider(X: float, eta: float) -> AxialFourthMixedJet:
            return fifth_provider(X, eta).fourth()

        def third_provider(X: float, eta: float) -> AxialThirdMixedJet:
            return fourth_provider(X, eta).third()

        return cls(
            order=order,
            phi_second_jet_provider=phi_second_provider,
            axial_third_mixed_jet_provider=third_provider,
            normalization_C=C,
            axial_fourth_mixed_jet_provider=fourth_provider,
            phi_third_mixed_jet_provider=phi_third_provider,
            axial_fifth_mixed_jet_provider=fifth_provider,
            phi_fourth_mixed_jet_provider=phi_fourth_provider,
            provenance=provenance
            or (
                "Lemma 5.2 compact repaired phi fourth-mixed/U fifth-mixed "
                "coefficient source; all lower angular/axial layers are exact "
                "projections and beta is derived only through analytic Eq. (5.2); "
                "formal-structure only."
            ),
        )


class Section5LowerHistoryPhiFourthMixedHierarchy(
    Section5LowerHistoryFifthMixedHierarchy
):
    """Lower history with fail-closed hierarchy-owned fourth-mixed phi data."""

    def _phi_fourth_value(
        self,
        order: int,
        X: float,
        eta: float,
    ) -> ProfileFourthMixedJet:
        source = self._source(order)
        if not isinstance(source, Section5PhiFourthMixedCoefficientJetSource):
            raise ValueError(
                f"coefficient order {source.order} does not own fourth-mixed phi data"
            )
        provider = source.phi_fourth_mixed_jet_provider
        if provider is None:  # defensive: __post_init__ already rejects this
            raise ValueError(
                f"coefficient order {source.order} does not own fourth-mixed phi data"
            )
        value = provider(float(X), float(eta))
        if not isinstance(value, ProfileFourthMixedJet):
            raise TypeError(
                "phi_fourth_mixed_jet_provider must return ProfileFourthMixedJet"
            )
        return value

    def phi_third_mixed_jet(
        self,
        order: int,
        X: float,
        eta: float,
    ) -> ProfileThirdMixedJet:
        """Return the fourth provider's exact third projection when available."""

        source = self._source(order)
        if not isinstance(source, Section5PhiFourthMixedCoefficientJetSource):
            return super().phi_third_mixed_jet(order, X, eta)

        fourth = self._phi_fourth_value(order, X, eta)
        projected = fourth.third()
        lower = super().phi_third_mixed_jet(order, X, eta)
        if projected != lower:
            raise ValueError(
                "fourth-mixed phi provider is incoherent with the hierarchy-owned "
                "third-mixed projection"
            )
        return projected

    def phi_fourth_mixed_jet(
        self,
        order: int,
        X: float,
        eta: float,
    ) -> ProfileFourthMixedJet:
        """Return the owned fourth-mixed phi jet and verify its lower projection."""

        value = self._phi_fourth_value(order, X, eta)
        if value.third() != self.phi_third_mixed_jet(order, X, eta):
            raise ValueError(
                "fourth-mixed phi provider is incoherent with the hierarchy-owned "
                "third-mixed projection"
            )
        return value
