"""Hierarchy-owned repaired sixth-mixed U and fifth-mixed beta jets for Section 5.

This layer closes the next ownership seam on the route to
``partial_eta^3(Omega/X)``. The coefficient source owns one authoritative
compact-repaired sixth-mixed ``U_n`` provider together with the already-landed
fifth-mixed ``phi_n`` provider. All lower U projections are checked for exact
coherence, and the fifth-mixed regular flux ``beta_n=V_n/X`` is derived only by
analytic Eq. (5.2); no caller-maintained beta table is accepted.

Genuine order-zero strong data and unrepaired base sixth-mixed jets remain Issue
#1/upstream inputs. This is fail-closed Stage-2 ``formal-structure``
infrastructure, not paper-exact coefficient materialization or convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .background_moment_repair_phi_fifth_mixed_jets import (
    BasePhiFifthMixedJetProvider,
    Lemma52RepairedPhiFifthMixedJetAdapter,
    ProfileFifthMixedJet,
)
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_moment_repair_sixth_mixed_jets import (
    BaseSixthMixedJetProvider,
    Lemma52RepairedSixthMixedJetAdapter,
)
from .background_preceding_diffusion_parameter_jet import ProfileThirdMixedJet
from .background_preceding_diffusion_second_parameter_jet import ProfileFourthMixedJet
from .background_regular_flux_fifth_mixed_jets import (
    AxialSixthMixedJet,
    RegularFluxFifthMixedJet,
    regular_flux_fifth_mixed_jet_eq_5_2,
)
from .background_regular_flux_fourth_mixed_jets import AxialFifthMixedJet
from .background_regular_flux_second_jets import AxialThirdMixedJet
from .background_regular_flux_third_mixed_jets import AxialFourthMixedJet
from .background_repaired_history_phi_fifth_mixed import (
    Section5LowerHistoryPhiFifthMixedHierarchy,
    Section5PhiFifthMixedCoefficientJetSource,
)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


@dataclass(frozen=True)
class Section5SixthMixedCoefficientJetSource(Section5PhiFifthMixedCoefficientJetSource):
    """Coefficient source with authoritative repaired sixth-mixed U data."""

    axial_sixth_mixed_jet_provider: BaseSixthMixedJetProvider | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if not callable(self.axial_sixth_mixed_jet_provider):
            raise TypeError("axial_sixth_mixed_jet_provider must be callable")

    @classmethod
    def from_lemma52_repair_full_sixth_mixed(
        cls, order: int, *, C: float, profile: Lemma52RepairedProfileAdapter,
        base_phi_fifth_mixed_jet: BasePhiFifthMixedJetProvider,
        base_U_sixth_mixed_jet: BaseSixthMixedJetProvider,
        provenance: str | None = None,
    ) -> "Section5SixthMixedCoefficientJetSource":
        """Own repaired phi-fifth and U-sixth jets for one positive order."""
        order = _positive_order(order)
        phi = Lemma52RepairedPhiFifthMixedJetAdapter(profile, C, base_phi_fifth_mixed_jet)
        axial = Lemma52RepairedSixthMixedJetAdapter(profile, base_U_sixth_mixed_jet)

        def phi_fifth_provider(X: float, eta: float) -> ProfileFifthMixedJet:
            return phi.phi_fifth_mixed_jet(X, eta)
        def phi_fourth_provider(X: float, eta: float) -> ProfileFourthMixedJet:
            return phi_fifth_provider(X, eta).fourth()
        def phi_third_provider(X: float, eta: float) -> ProfileThirdMixedJet:
            return phi_fourth_provider(X, eta).third()
        def phi_second_provider(X: float, eta: float):
            return phi_third_provider(X, eta).second()
        def sixth_provider(X: float, eta: float) -> AxialSixthMixedJet:
            return axial.U_sixth_mixed_jet(X, eta)
        def fifth_provider(X: float, eta: float) -> AxialFifthMixedJet:
            return sixth_provider(X, eta).fifth()
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
            phi_fifth_mixed_jet_provider=phi_fifth_provider,
            axial_sixth_mixed_jet_provider=sixth_provider,
            provenance=provenance or (
                "Lemma 5.2 compact repaired phi fifth-mixed/U sixth-mixed coefficient "
                "source; all lower layers are exact projections and beta is derived "
                "only through analytic Eq. (5.2); formal-structure only."
            ),
        )


class Section5LowerHistorySixthMixedHierarchy(Section5LowerHistoryPhiFifthMixedHierarchy):
    """Lower history with fail-closed hierarchy-owned sixth-mixed U data."""

    def _sixth_value(self, order: int, X: float, eta: float) -> AxialSixthMixedJet:
        source = self._source(order)
        if not isinstance(source, Section5SixthMixedCoefficientJetSource):
            raise ValueError(f"coefficient order {source.order} does not own sixth-mixed U data")
        provider = source.axial_sixth_mixed_jet_provider
        if provider is None:
            raise ValueError(f"coefficient order {source.order} does not own sixth-mixed U data")
        value = provider(float(X), float(eta))
        if not isinstance(value, AxialSixthMixedJet):
            raise TypeError("axial_sixth_mixed_jet_provider must return AxialSixthMixedJet")
        return value

    def axial_fifth_mixed_jet(self, order: int, X: float, eta: float) -> AxialFifthMixedJet:
        """Return the sixth provider's exact fifth projection when available."""
        source = self._source(order)
        if not isinstance(source, Section5SixthMixedCoefficientJetSource):
            return super().axial_fifth_mixed_jet(order, X, eta)
        sixth = self._sixth_value(order, X, eta)
        projected = sixth.fifth()
        lower = super().axial_fifth_mixed_jet(order, X, eta)
        if projected != lower:
            raise ValueError(
                "sixth-mixed U provider is incoherent with the hierarchy-owned fifth-mixed projection"
            )
        return projected

    def axial_sixth_mixed_jet(self, order: int, X: float, eta: float) -> AxialSixthMixedJet:
        """Return owned sixth-mixed U and verify its lower projection."""
        value = self._sixth_value(order, X, eta)
        if value.fifth() != self.axial_fifth_mixed_jet(order, X, eta):
            raise ValueError(
                "sixth-mixed U provider is incoherent with the hierarchy-owned fifth-mixed projection"
            )
        return value

    def beta_fifth_mixed_jet(self, order: int, X: float, eta: float) -> RegularFluxFifthMixedJet:
        """Derive owned fifth-mixed ``beta=V/X`` through analytic Eq. (5.2)."""
        self._source(order)
        def axial_provider(x: float, e: float) -> AxialSixthMixedJet:
            return self.axial_sixth_mixed_jet(order, x, e)
        return regular_flux_fifth_mixed_jet_eq_5_2(
            self.h, order, X, eta, axial_provider, quadrature_points=self.quadrature_points,
        )
