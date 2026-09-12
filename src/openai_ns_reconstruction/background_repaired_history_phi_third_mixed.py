"""Hierarchy-owned third-mixed ``phi_n`` jets for Section 5 lower history.

The base :mod:`background_repaired_history` hierarchy already owns the repaired
``phi_n`` second jet and, on its stronger path, the repaired fourth-mixed
``U_n`` jet together with the Eq. (5.2) third-mixed regular flux.  The analytic
eta derivative of the strict-lower ``precedingDiffusion`` term needs one more
layer from ``phi_n``: ``phi_XXeta``, ``phi_Xetaeta`` and ``phi_etaetaeta``.

This module adds that layer without changing the compatibility surface of the
older hierarchy.  A strong coefficient source stores one actual
``Lemma52RepairedPhiThirdMixedJetAdapter`` provider and projects its first six
fields to the inherited second-jet provider.  The matching hierarchy subclass
therefore cannot pair a third-mixed phi record with an unrelated lower-order
phi record: every query checks that the projection is exactly the ordinary
hierarchy-owned second jet and fails closed otherwise.

The Lemma 5.2 constructor is positive-order only.  Order zero may still be
materialized explicitly with this source type once Issue #1 supplies a genuine
leading third-mixed profile.  Until then, missing order-zero third-mixed data are
a deliberate blocker for a fully hierarchy-owned ``partial_eta f_1`` path.

Upstream unrepaired jets, moment/patch eta-jets, normalization C and the leading
profile remain inputs.  This is Stage-2 ``formal-structure`` infrastructure and
must not be interpreted as a paper-exact recursive coefficient or velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .background_moment_repair_fourth_mixed_jets import (
    BaseFourthMixedJetProvider,
    Lemma52RepairedFourthMixedJetAdapter,
)
from .background_moment_repair_phi_third_mixed_jets import (
    BasePhiThirdMixedJetProvider,
    Lemma52RepairedPhiThirdMixedJetAdapter,
)
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_preceding_diffusion_parameter_jet import ProfileThirdMixedJet
from .background_regular_flux_second_jets import AxialThirdMixedJet
from .background_regular_flux_third_mixed_jets import AxialFourthMixedJet
from .background_repaired_history import (
    Section5CoefficientJetSource,
    Section5LowerHistoryJetHierarchy,
)


def _positive_order(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


@dataclass(frozen=True)
class Section5PhiThirdMixedCoefficientJetSource(Section5CoefficientJetSource):
    """Coefficient source whose third-mixed phi jet is hierarchy-owned.

    ``phi_third_mixed_jet_provider`` is required.  The Lemma-5.2 factory below
    also makes the fourth-mixed U provider authoritative, so the strong phi and
    strong U/beta paths live in one coefficient source rather than in parallel
    caller-maintained tables.
    """

    phi_third_mixed_jet_provider: BasePhiThirdMixedJetProvider | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if not callable(self.phi_third_mixed_jet_provider):
            raise TypeError("phi_third_mixed_jet_provider must be callable")

    @classmethod
    def from_lemma52_repair_full_mixed(
        cls,
        order: int,
        *,
        C: float,
        profile: Lemma52RepairedProfileAdapter,
        base_phi_third_mixed_jet: BasePhiThirdMixedJetProvider,
        base_U_fourth_mixed_jet: BaseFourthMixedJetProvider,
        provenance: str | None = None,
    ) -> "Section5PhiThirdMixedCoefficientJetSource":
        """Own repaired phi third-mixed and U fourth-mixed jets for one order.

        The same repaired phi provider supplies both derivative layers: the
        inherited second-jet path is its exact projection.  Likewise the same
        repaired fourth-mixed U provider supplies the inherited third-mixed U
        path.  No independent beta derivative provider is accepted; the base
        hierarchy continues to derive beta only through analytic Eq. (5.2).
        """

        order = _positive_order(order)
        phi = Lemma52RepairedPhiThirdMixedJetAdapter(
            profile,
            C,
            base_phi_third_mixed_jet,
        )
        axial = Lemma52RepairedFourthMixedJetAdapter(
            profile,
            base_U_fourth_mixed_jet,
        )

        def phi_third_provider(X: float, eta: float) -> ProfileThirdMixedJet:
            return phi.phi_third_mixed_jet(X, eta)

        def phi_second_provider(X: float, eta: float):
            return phi_third_provider(X, eta).second()

        def fourth_provider(X: float, eta: float) -> AxialFourthMixedJet:
            return axial.U_fourth_mixed_jet(X, eta)

        def third_provider(X: float, eta: float) -> AxialThirdMixedJet:
            return fourth_provider(X, eta).third()

        return cls(
            order=order,
            phi_second_jet_provider=phi_second_provider,
            axial_third_mixed_jet_provider=third_provider,
            normalization_C=C,
            axial_fourth_mixed_jet_provider=fourth_provider,
            phi_third_mixed_jet_provider=phi_third_provider,
            provenance=provenance
            or (
                "Lemma 5.2 compact repaired phi third-mixed/U fourth-mixed "
                "coefficient source; lower phi/U layers are exact projections "
                "and beta is derived only through analytic Eq. (5.2); "
                "formal-structure only."
            ),
        )


class Section5LowerHistoryPhiThirdMixedHierarchy(Section5LowerHistoryJetHierarchy):
    """Lower-history hierarchy with a fail-closed owned phi third-mixed path."""

    def phi_third_mixed_jet(
        self,
        order: int,
        X: float,
        eta: float,
    ) -> ProfileThirdMixedJet:
        """Return the owned strong phi jet and verify its lower projection.

        Older coefficient sources intentionally fail here.  This makes the
        absence of Issue-#1 leading third-mixed data explicit instead of
        silently substituting a caller-supplied derivative table.
        """

        source = self._source(order)
        if not isinstance(source, Section5PhiThirdMixedCoefficientJetSource):
            raise ValueError(
                f"coefficient order {source.order} does not own third-mixed phi data"
            )
        provider = source.phi_third_mixed_jet_provider
        if provider is None:  # defensive: __post_init__ already rejects this
            raise ValueError(
                f"coefficient order {source.order} does not own third-mixed phi data"
            )
        value = provider(float(X), float(eta))
        if not isinstance(value, ProfileThirdMixedJet):
            raise TypeError(
                "phi_third_mixed_jet_provider must return ProfileThirdMixedJet"
            )
        if value.second() != self.phi_second_jet(order, X, eta):
            raise ValueError(
                "third-mixed phi provider is incoherent with the hierarchy-owned "
                "second-jet projection"
            )
        return value
