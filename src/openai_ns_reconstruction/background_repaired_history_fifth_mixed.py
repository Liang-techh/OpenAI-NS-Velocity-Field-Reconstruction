"""Hierarchy-owned repaired fifth-mixed U jets for Section 5.

The landed Lemma 5.2 adapter already produces the three fifth-mixed U
derivatives needed by the analytic Eq. (5.2) fourth-mixed regular-flux path.
This module makes that stronger U layer part of the same coefficient source
that already owns repaired phi third-mixed data.  The fifth-mixed provider is
authoritative: its fourth- and third-mixed projections are the inherited U
providers, and every hierarchy query checks exact projection coherence.

No independent beta table is accepted.  ``beta_fourth_mixed_jet`` is always
derived from the owned fifth-mixed U provider through the landed analytic
Eq. (5.2) adapter.  Genuine order-zero strong data remain an Issue #1 input, so
this is fail-closed Stage-2 ``formal-structure`` infrastructure rather than a
paper-exact recursive coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .background_moment_repair_fifth_mixed_jets import (
    BaseFifthMixedJetProvider,
    Lemma52RepairedFifthMixedJetAdapter,
)
from .background_moment_repair_phi_third_mixed_jets import (
    BasePhiThirdMixedJetProvider,
    Lemma52RepairedPhiThirdMixedJetAdapter,
)
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_omega_second_parameter_jet import RegularFluxFourthMixedJet
from .background_preceding_diffusion_parameter_jet import ProfileThirdMixedJet
from .background_regular_flux_fourth_mixed_jets import (
    AxialFifthMixedJet,
    regular_flux_fourth_mixed_jet_eq_5_2,
)
from .background_regular_flux_second_jets import AxialThirdMixedJet
from .background_regular_flux_third_mixed_jets import AxialFourthMixedJet
from .background_repaired_history_phi_third_mixed import (
    Section5LowerHistoryPhiThirdMixedHierarchy,
    Section5PhiThirdMixedCoefficientJetSource,
)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


@dataclass(frozen=True)
class Section5FifthMixedCoefficientJetSource(
    Section5PhiThirdMixedCoefficientJetSource
):
    """Coefficient source with one authoritative repaired fifth-mixed U provider."""

    axial_fifth_mixed_jet_provider: BaseFifthMixedJetProvider | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if not callable(self.axial_fifth_mixed_jet_provider):
            raise TypeError("axial_fifth_mixed_jet_provider must be callable")

    @classmethod
    def from_lemma52_repair_full_fifth_mixed(
        cls,
        order: int,
        *,
        C: float,
        profile: Lemma52RepairedProfileAdapter,
        base_phi_third_mixed_jet: BasePhiThirdMixedJetProvider,
        base_U_fifth_mixed_jet: BaseFifthMixedJetProvider,
        provenance: str | None = None,
    ) -> "Section5FifthMixedCoefficientJetSource":
        """Own repaired phi-third and U-fifth jets for one positive order.

        The repaired U fifth-mixed adapter is the only strong U source.
        Fourth- and third-mixed providers are exact projections of it, while
        beta remains derived only through the analytic Eq. (5.2) hierarchy
        methods.  Lemma 5.2 is positive-order only, so order zero must still be
        supplied explicitly by Issue #1.
        """

        order = _positive_order(order)
        phi = Lemma52RepairedPhiThirdMixedJetAdapter(
            profile,
            C,
            base_phi_third_mixed_jet,
        )
        axial = Lemma52RepairedFifthMixedJetAdapter(
            profile,
            base_U_fifth_mixed_jet,
        )

        def phi_third_provider(X: float, eta: float) -> ProfileThirdMixedJet:
            return phi.phi_third_mixed_jet(X, eta)

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
            provenance=provenance
            or (
                "Lemma 5.2 compact repaired phi third-mixed/U fifth-mixed "
                "coefficient source; lower U layers are exact projections and "
                "beta is derived only through analytic Eq. (5.2); "
                "formal-structure only."
            ),
        )


class Section5LowerHistoryFifthMixedHierarchy(
    Section5LowerHistoryPhiThirdMixedHierarchy
):
    """Lower history with fail-closed hierarchy-owned fifth-mixed U data."""

    def _fifth_value(
        self,
        order: int,
        X: float,
        eta: float,
    ) -> AxialFifthMixedJet:
        source = self._source(order)
        if not isinstance(source, Section5FifthMixedCoefficientJetSource):
            raise ValueError(
                f"coefficient order {source.order} does not own fifth-mixed U data"
            )
        provider = source.axial_fifth_mixed_jet_provider
        if provider is None:  # defensive: __post_init__ already rejects this
            raise ValueError(
                f"coefficient order {source.order} does not own fifth-mixed U data"
            )
        value = provider(float(X), float(eta))
        if not isinstance(value, AxialFifthMixedJet):
            raise TypeError(
                "axial_fifth_mixed_jet_provider must return AxialFifthMixedJet"
            )
        return value

    def axial_fourth_mixed_jet(
        self,
        order: int,
        X: float,
        eta: float,
    ) -> AxialFourthMixedJet:
        """Return the fifth provider's exact fourth projection when available."""

        source = self._source(order)
        if not isinstance(source, Section5FifthMixedCoefficientJetSource):
            return super().axial_fourth_mixed_jet(order, X, eta)

        fifth = self._fifth_value(order, X, eta)
        projected = fifth.fourth()
        provider = source.axial_fourth_mixed_jet_provider
        if provider is None:
            raise ValueError(
                f"coefficient order {source.order} does not own fourth-mixed U data"
            )
        lower = provider(float(X), float(eta))
        if not isinstance(lower, AxialFourthMixedJet):
            raise TypeError(
                "axial_fourth_mixed_jet_provider must return AxialFourthMixedJet"
            )
        if projected != lower:
            raise ValueError(
                "fifth-mixed U provider is incoherent with the hierarchy-owned "
                "fourth-mixed projection"
            )
        return projected

    def axial_fifth_mixed_jet(
        self,
        order: int,
        X: float,
        eta: float,
    ) -> AxialFifthMixedJet:
        """Return the owned fifth-mixed U jet and verify its lower projection."""

        value = self._fifth_value(order, X, eta)
        if value.fourth() != self.axial_fourth_mixed_jet(order, X, eta):
            raise ValueError(
                "fifth-mixed U provider is incoherent with the hierarchy-owned "
                "fourth-mixed projection"
            )
        return value

    def beta_fourth_mixed_jet(
        self,
        order: int,
        X: float,
        eta: float,
    ) -> RegularFluxFourthMixedJet:
        """Derive the owned fourth-mixed ``beta=V/X`` jet via Eq. (5.2)."""

        self._source(order)

        def axial_provider(x: float, e: float) -> AxialFifthMixedJet:
            return self.axial_fifth_mixed_jet(order, x, e)

        return regular_flux_fourth_mixed_jet_eq_5_2(
            self.h,
            order,
            X,
            eta,
            axial_provider,
            quadrature_points=self.quadrature_points,
        )
