"""Hierarchy-owned lower-history jets for the Section 5 PositiveAxis recursion.

The landed Lemma 5.2 adapters already provide three pieces separately:

* the repaired PositiveAxis angular variable ``phi_n`` through second
  ``(X, eta)`` order;
* the repaired axial coefficient ``U_n`` through the mixed third derivatives
  required by Eq. (5.2); and
* the analytic Eq. (5.2) map from that axial jet to the regular flux
  ``beta_n = V_n / X`` through second order.

A later derivative layer also provides the repaired axial fourth-mixed jet and
the corresponding third-mixed ``beta_n`` jet required by the analytic
``partial_eta(Omega/X)`` path.  This module makes those pieces hierarchy-owned
when the stronger axial jet is available, while keeping the older third-mixed
source path valid for coefficients whose upstream data stop there.

The strict-lower-history solver consumes order-indexed provider functions.  This
module makes one coefficient hierarchy own those providers and connects them
directly to the landed ``actualLowerSource`` / Eq. (5.7) bridge.  It deliberately
keeps order zero generic because the leading profile belongs to Issue #1, while
positive orders may be materialized with the actual Lemma 5.2 compact repair.

No coefficient is solved here.  Upstream base jets, moment/patch eta-jets, the
normalization C, and the leading order remain explicit data.  Consequently this
is solver infrastructure / ``formal-structure`` only; it never promotes the
background hierarchy or final velocity to paper-exact status.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Callable, Sequence
import math

import numpy as np

from .background_lower_history_solver import (
    LowerHistoryPositiveAxisProviders,
    first_positive_order_term_from_lower_history,
    lower_history_positive_axis_providers,
    positive_axis_eq_5_7_fields_from_lower_history,
)
from .background_lower_history_source import ProfileSecondJet
from .background_moment_repair_fourth_mixed_jets import (
    BaseFourthMixedJetProvider,
    Lemma52RepairedFourthMixedJetAdapter,
)
from .background_moment_repair_phi_jets import (
    BasePhiSecondJetProvider,
    Lemma52RepairedPhiSecondJetAdapter,
)
from .background_moment_repair_profile import Lemma52RepairedProfileAdapter
from .background_moment_repair_third_mixed_jets import (
    BaseThirdMixedJetProvider,
    Lemma52RepairedThirdMixedJetAdapter,
)
from .background_omega_parameter_jet import RegularFluxThirdMixedJet
from .background_positive_axis import PositiveAxisFields
from .background_regular_flux_second_jets import (
    AxialThirdMixedJet,
    AxialThirdMixedJetProvider,
    regular_flux_second_jet_eq_5_2,
)
from .background_regular_flux_third_mixed_jets import (
    AxialFourthMixedJet,
    AxialFourthMixedJetProvider,
    regular_flux_third_mixed_jet_eq_5_2,
)
from .coordinates import validate_h


SecondJetProvider = Callable[[float, float], ProfileSecondJet]


def _order(value: int, *, positive: bool = False) -> int:
    lower = 1 if positive else 0
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < lower:
        qualifier = "positive " if positive else "nonnegative "
        raise ValueError(f"order must be a {qualifier}integer")
    return int(value)


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _second_from_third(jet: AxialThirdMixedJet) -> ProfileSecondJet:
    if not isinstance(jet, AxialThirdMixedJet):
        raise TypeError("axial_third_mixed_jet_provider must return AxialThirdMixedJet")
    return ProfileSecondJet(
        value=jet.value,
        radial=jet.radial,
        radial2=jet.radial2,
        parameter=jet.parameter,
        radial_parameter=jet.radial_parameter,
        parameter2=jet.parameter2,
    )


@dataclass(frozen=True)
class Section5CoefficientJetSource:
    """One hierarchy-owned coefficient jet source.

    ``phi_second_jet_provider`` supplies the PositiveAxis ``phi_n`` jet.
    ``axial_third_mixed_jet_provider`` supplies a coherent repaired/unrepaired
    ``U_n`` jet and is also the source used to construct ``beta_n`` through
    second order via Eq. (5.2).

    When ``axial_fourth_mixed_jet_provider`` is present it is authoritative for
    the axial coefficient: the hierarchy projects its first nine fields back to
    the third-mixed jet and derives the third-mixed regular-flux jet only through
    the analytic Eq. (5.2) adapter.  This prevents the stronger derivative layer
    from being paired with an unrelated lower-order U record.

    ``normalization_C`` is populated by the Lemma-5.2 constructors so a
    hierarchy can reject accidental mixing of repairs built with a different C.
    Generic sources keep it ``None``; this is required for the Issue-#1 leading
    coefficient until that upstream materialization is complete.
    """

    order: int
    phi_second_jet_provider: SecondJetProvider
    axial_third_mixed_jet_provider: AxialThirdMixedJetProvider
    provenance: str
    normalization_C: float | None = None
    axial_fourth_mixed_jet_provider: AxialFourthMixedJetProvider | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "order", _order(self.order))
        if not callable(self.phi_second_jet_provider):
            raise TypeError("phi_second_jet_provider must be callable")
        if not callable(self.axial_third_mixed_jet_provider):
            raise TypeError("axial_third_mixed_jet_provider must be callable")
        if self.axial_fourth_mixed_jet_provider is not None and not callable(
            self.axial_fourth_mixed_jet_provider
        ):
            raise TypeError("axial_fourth_mixed_jet_provider must be callable")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")
        if self.normalization_C is not None:
            object.__setattr__(
                self,
                "normalization_C",
                _positive(self.normalization_C, "normalization_C"),
            )

    @classmethod
    def from_lemma52_repair(
        cls,
        order: int,
        *,
        C: float,
        profile: Lemma52RepairedProfileAdapter,
        base_phi_second_jet: BasePhiSecondJetProvider,
        base_U_third_mixed_jet: BaseThirdMixedJetProvider,
        provenance: str | None = None,
    ) -> "Section5CoefficientJetSource":
        """Own the actual compact Lemma-5.2 repaired jets for one positive order.

        This compatibility constructor materializes the derivative layer needed
        by the ordinary second-jet lower-history solve.  Use
        :meth:`from_lemma52_repair_fourth_mixed` when the upstream coefficient
        also supplies the fourth-mixed U data needed by the analytic
        ``partial_eta(Omega/X)`` route.

        Lemma 5.2 belongs to the positive-order Section-5 correction hierarchy,
        so order zero is rejected rather than being silently treated as a
        repaired leading profile.
        """

        order = _order(order, positive=True)
        C = _positive(C, "C")
        if not isinstance(profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        phi = Lemma52RepairedPhiSecondJetAdapter(profile, C, base_phi_second_jet)
        axial = Lemma52RepairedThirdMixedJetAdapter(profile, base_U_third_mixed_jet)
        return cls(
            order=order,
            phi_second_jet_provider=phi.phi_second_jet,
            axial_third_mixed_jet_provider=axial.U_third_mixed_jet,
            normalization_C=C,
            provenance=provenance
            or (
                "Lemma 5.2 compact repaired phi/U coefficient source; beta is "
                "derived only through analytic Eq. (5.2); formal-structure only."
            ),
        )

    @classmethod
    def from_lemma52_repair_fourth_mixed(
        cls,
        order: int,
        *,
        C: float,
        profile: Lemma52RepairedProfileAdapter,
        base_phi_second_jet: BasePhiSecondJetProvider,
        base_U_fourth_mixed_jet: BaseFourthMixedJetProvider,
        provenance: str | None = None,
    ) -> "Section5CoefficientJetSource":
        """Own one repaired coefficient through the fourth-mixed axial layer.

        The fourth-mixed adapter is the single authoritative U provider.  Its
        projection supplies the already-landed third-mixed path, while the full
        jet supplies ``beta_XXeta``, ``beta_Xetaeta`` and ``beta_etaetaeta`` via
        the analytic Eq. (5.2) adapter.  No independently supplied beta jet or
        duplicate third-mixed U table is accepted by this constructor.
        """

        order = _order(order, positive=True)
        C = _positive(C, "C")
        if not isinstance(profile, Lemma52RepairedProfileAdapter):
            raise TypeError("profile must be a Lemma52RepairedProfileAdapter")
        if not callable(base_U_fourth_mixed_jet):
            raise TypeError("base_U_fourth_mixed_jet must be callable")

        phi = Lemma52RepairedPhiSecondJetAdapter(profile, C, base_phi_second_jet)
        axial = Lemma52RepairedFourthMixedJetAdapter(
            profile, base_U_fourth_mixed_jet
        )

        def fourth_provider(X: float, eta: float) -> AxialFourthMixedJet:
            return axial.U_fourth_mixed_jet(X, eta)

        def third_provider(X: float, eta: float) -> AxialThirdMixedJet:
            return fourth_provider(X, eta).third()

        return cls(
            order=order,
            phi_second_jet_provider=phi.phi_second_jet,
            axial_third_mixed_jet_provider=third_provider,
            normalization_C=C,
            axial_fourth_mixed_jet_provider=fourth_provider,
            provenance=provenance
            or (
                "Lemma 5.2 compact repaired phi/U coefficient source through "
                "fourth-mixed U; beta second/third-mixed jets are derived only "
                "through analytic Eq. (5.2); formal-structure only."
            ),
        )


@dataclass(frozen=True)
class Section5LowerHistoryJetHierarchy:
    """Contiguous coefficient history consumed by ``actualLowerSource``.

    Sources must be exactly orders ``0,...,N``.  A request for recursive order
    ``n`` is admitted only when the hierarchy owns every strict lower order
    ``0,...,n-1``.  The regular flux is generated on demand from each owned
    axial jet with Eq. (5.2); no independent beta provider is accepted.

    The ordinary lower-history path needs only third-mixed U data.  The stronger
    ``beta_third_mixed_jet`` path fails closed unless that coefficient source
    owns a fourth-mixed U provider.
    """

    h: float
    C: float
    sources: Sequence[Section5CoefficientJetSource]
    quadrature_points: int = 32
    provenance: str = (
        "Contiguous Section-5 lower-history provider assembled from owned phi/U "
        "jets with beta derived by Eq. (5.2); formal-structure only."
    )

    def __post_init__(self) -> None:
        object.__setattr__(self, "h", validate_h(self.h))
        object.__setattr__(self, "C", _positive(self.C, "C"))
        if (
            isinstance(self.quadrature_points, bool)
            or not isinstance(self.quadrature_points, Integral)
            or int(self.quadrature_points) < 1
        ):
            raise ValueError("quadrature_points must be a positive integer")
        object.__setattr__(self, "quadrature_points", int(self.quadrature_points))
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty string")

        sources = tuple(self.sources)
        if not sources:
            raise ValueError("sources must contain at least the order-zero coefficient")
        if not all(isinstance(source, Section5CoefficientJetSource) for source in sources):
            raise TypeError("sources must contain Section5CoefficientJetSource entries")
        orders = tuple(source.order for source in sources)
        expected = tuple(range(len(sources)))
        if orders != expected:
            raise ValueError(
                f"sources must be contiguous and ordered as {expected}; got {orders}"
            )
        for source in sources:
            if source.normalization_C is not None and not math.isclose(
                source.normalization_C, self.C, rel_tol=0.0, abs_tol=0.0
            ):
                raise ValueError(
                    f"source order {source.order} was repaired with a different C"
                )
        object.__setattr__(self, "sources", sources)

    def _source(self, order: int) -> Section5CoefficientJetSource:
        order = _order(order)
        if order >= len(self.sources):
            raise ValueError(f"coefficient order {order} is not owned by this hierarchy")
        return self.sources[order]

    def phi_second_jet(self, order: int, X: float, eta: float) -> ProfileSecondJet:
        """Return the owned ``phi_order`` second jet."""

        value = self._source(order).phi_second_jet_provider(float(X), float(eta))
        if not isinstance(value, ProfileSecondJet):
            raise TypeError("phi_second_jet_provider must return ProfileSecondJet")
        return value

    def axial_fourth_mixed_jet(
        self, order: int, X: float, eta: float
    ) -> AxialFourthMixedJet:
        """Return the owned fourth-mixed axial jet, or fail closed if absent."""

        source = self._source(order)
        provider = source.axial_fourth_mixed_jet_provider
        if provider is None:
            raise ValueError(
                f"coefficient order {source.order} does not own fourth-mixed U data"
            )
        value = provider(float(X), float(eta))
        if not isinstance(value, AxialFourthMixedJet):
            raise TypeError(
                "axial_fourth_mixed_jet_provider must return AxialFourthMixedJet"
            )
        return value

    def axial_third_mixed_jet(
        self, order: int, X: float, eta: float
    ) -> AxialThirdMixedJet:
        """Return the owned coherent axial jet used by both U and beta paths.

        If a fourth-mixed provider exists it is authoritative and projected to
        third order here, so the lower-history path cannot silently diverge from
        the stronger derivative layer.
        """

        source = self._source(order)
        if source.axial_fourth_mixed_jet_provider is not None:
            return self.axial_fourth_mixed_jet(order, X, eta).third()
        value = source.axial_third_mixed_jet_provider(float(X), float(eta))
        if not isinstance(value, AxialThirdMixedJet):
            raise TypeError(
                "axial_third_mixed_jet_provider must return AxialThirdMixedJet"
            )
        return value

    def axial_second_jet(self, order: int, X: float, eta: float) -> ProfileSecondJet:
        """Project the owned coherent axial third-mixed jet to second order."""

        return _second_from_third(self.axial_third_mixed_jet(order, X, eta))

    def beta_second_jet(self, order: int, X: float, eta: float) -> ProfileSecondJet:
        """Derive the owned ``beta_order=V_order/X`` second jet via Eq. (5.2)."""

        self._source(order)

        def axial_provider(x: float, e: float) -> AxialThirdMixedJet:
            return self.axial_third_mixed_jet(order, x, e)

        return regular_flux_second_jet_eq_5_2(
            self.h,
            order,
            X,
            eta,
            axial_provider,
            quadrature_points=self.quadrature_points,
        )

    def beta_third_mixed_jet(
        self, order: int, X: float, eta: float
    ) -> RegularFluxThirdMixedJet:
        """Derive the owned third-mixed ``beta=V/X`` jet via Eq. (5.2).

        This path intentionally has no fallback to caller-supplied beta
        derivatives.  It requires the same coefficient source to own the
        fourth-mixed U jet and fails closed otherwise.
        """

        self._source(order)

        def axial_provider(x: float, e: float) -> AxialFourthMixedJet:
            return self.axial_fourth_mixed_jet(order, x, e)

        return regular_flux_third_mixed_jet_eq_5_2(
            self.h,
            order,
            X,
            eta,
            axial_provider,
            quadrature_points=self.quadrature_points,
        )

    def positive_axis_providers(self, order: int) -> LowerHistoryPositiveAxisProviders:
        """Connect owned strict lower history directly to ``actualLowerSource``.

        Only orders ``0,...,order-1`` are queried by the landed bridge.  The
        current order need not exist yet, which is exactly the recursive use
        case: solve coefficient n from already materialized lower coefficients.
        """

        order = _order(order, positive=True)
        if order > len(self.sources):
            raise ValueError(
                "requested recursive order is missing one or more strict lower coefficients"
            )
        return lower_history_positive_axis_providers(
            self.h,
            order,
            self.phi_second_jet,
            self.axial_second_jet,
            self.beta_second_jet,
        )

    def positive_axis_fields(self, order: int) -> PositiveAxisFields:
        """Return the genuine Eq. (5.7) fields generated from owned lower history."""

        order = _order(order, positive=True)
        if order > len(self.sources):
            raise ValueError(
                "requested recursive order is missing one or more strict lower coefficients"
            )
        return positive_axis_eq_5_7_fields_from_lower_history(
            self.h,
            order,
            self.C,
            self.phi_second_jet,
            self.axial_second_jet,
            self.beta_second_jet,
        )

    def first_picard_term(
        self,
        order: int,
        xi: float,
        eta: float,
        *,
        quadrature_points: int | None = None,
    ) -> np.ndarray:
        """Evaluate ``G f_n`` using only hierarchy-owned strict lower data."""

        order = _order(order, positive=True)
        if order > len(self.sources):
            raise ValueError(
                "requested recursive order is missing one or more strict lower coefficients"
            )
        q = self.quadrature_points if quadrature_points is None else quadrature_points
        if isinstance(q, bool) or not isinstance(q, Integral) or int(q) < 1:
            raise ValueError("quadrature_points must be a positive integer")
        return first_positive_order_term_from_lower_history(
            order,
            xi,
            eta,
            h=self.h,
            C=self.C,
            phi_provider=self.phi_second_jet,
            axial_provider=self.axial_second_jet,
            beta_provider=self.beta_second_jet,
            quadrature_points=int(q),
        )
