"""Finite Picard family bridge for the formal natural-axis coefficient graph.

The bridge starts from the actual reference pair ``x0`` and applies
``T(x) = x0 + R(x)/(2*Lambda)`` a requested finite number of times.  Each
application uses the same remainder graph as
:class:`FormalAxisCoefficientSolverState`; the solver's optional input-family
path keeps nonlinear operands distinct from the updated output family.

This module exposes only finite angular/axial coefficient prefixes.  It does
not materialize pressure for an updated iterate, establish a fixed point, or
provide a global weighted-space certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from functools import lru_cache
import math

from .axis_coefficient_formal_solver import FormalAxisCoefficientSolverState
from .axis_coefficient_mixed_scale import MixedScaleCoefficient, MixedScaleFamily


_DECIMAL_PRECISION = 96
_WINDOW_LEFT = -11.0 / 10.0
_WINDOW_RIGHT = 11.0 / 10.0


def _index(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _eta_in_window(value: float) -> float:
    eta = float(value)
    if not math.isfinite(eta) or not _WINDOW_LEFT <= eta <= _WINDOW_RIGHT:
        raise ValueError("eta must be finite and lie in the pinned window [-11/10,11/10]")
    return eta


def _reference_pair_families(
    solver: FormalAxisCoefficientSolverState,
    eta: float,
) -> tuple[MixedScaleFamily, MixedScaleFamily]:
    """Build one local cached pair of Decimal reference families."""

    eta = _eta_in_window(eta)

    @lru_cache(maxsize=None)
    def reference_pair(n: int, m: int) -> tuple[Decimal, Decimal]:
        return solver.reference_jet_pair(n, m, eta)

    @lru_cache(maxsize=None)
    def reference_phi(n: int, m: int) -> MixedScaleCoefficient:
        phi, _ = reference_pair(n, m)
        return MixedScaleCoefficient.channel(0, 0, phi)

    @lru_cache(maxsize=None)
    def reference_u(n: int, m: int) -> MixedScaleCoefficient:
        _, u = reference_pair(n, m)
        return MixedScaleCoefficient.channel(0, 0, u)

    def _local_eta(z: float) -> None:
        if _eta_in_window(z) != eta:
            raise ValueError("finite Picard family graph requires one fixed eta")

    def phi_family(n: int, m: int, z: float) -> MixedScaleCoefficient:
        _local_eta(z)
        return reference_phi(_index(n, "n"), _index(m, "m"))

    def u_family(n: int, m: int, z: float) -> MixedScaleCoefficient:
        _local_eta(z)
        return reference_u(_index(n, "n"), _index(m, "m"))

    return phi_family, u_family


@dataclass(frozen=True)
class FormalAxisPicardFamilyState:
    """Actual-solver-anchored finite Picard family evaluator."""

    solver: FormalAxisCoefficientSolverState

    def __post_init__(self) -> None:
        if not isinstance(self.solver, FormalAxisCoefficientSolverState):
            raise TypeError("solver must be FormalAxisCoefficientSolverState")

    @property
    def Lambda(self) -> Decimal:
        return self.solver.Lambda

    @property
    def epsilon(self) -> float:
        return self.solver.epsilon

    @property
    def finite_picard_materialized(self) -> bool:
        return True

    @property
    def formal_coefficients_materialized(self) -> bool:
        return True

    @property
    def fixed_point_materialized(self) -> bool:
        return False

    @property
    def fixed_point_convergence_certified(self) -> bool:
        return False

    @property
    def global_axis_norm_certified(self) -> bool:
        return False

    @property
    def truncation_certified(self) -> bool:
        return False

    @property
    def paper_exact(self) -> bool:
        return False

    def jet_prefix(
        self,
        iterations: int,
        max_n: int,
        m: int,
        eta: float,
    ) -> tuple[tuple[MixedScaleCoefficient, MixedScaleCoefficient], ...]:
        """Return rows ``0..max_n`` after ``iterations`` applications of ``T``.

        ``iterations = 0`` is the actual reference pair ``x0``.  Every map
        application builds one local cached graph, and no graph cache is kept
        on this state between calls.
        """

        iterations = _index(iterations, "iterations")
        max_n = _index(max_n, "max_n")
        m = _index(m, "m")
        eta = _eta_in_window(eta)

        with localcontext() as context:
            context.prec = _DECIMAL_PRECISION
            phi_family, u_family = _reference_pair_families(self.solver, eta)
            for _ in range(iterations):
                phi_family, u_family = self.solver.picard_map_families(
                    (phi_family, u_family),
                    eta,
                )
            return tuple(
                (phi_family(n, m, eta), u_family(n, m, eta))
                for n in range(max_n + 1)
            )

    def jet_pair(
        self,
        iterations: int,
        n: int,
        m: int,
        eta: float,
    ) -> tuple[MixedScaleCoefficient, MixedScaleCoefficient]:
        """Return one finite-iterate row from :meth:`jet_prefix`."""

        n = _index(n, "n")
        return self.jet_prefix(iterations, n, m, eta)[n]


def formal_axis_picard_family_state(
    solver: FormalAxisCoefficientSolverState,
) -> FormalAxisPicardFamilyState:
    """Bind a finite Picard family bridge to one genuine solver anchor."""

    return FormalAxisPicardFamilyState(solver=solver)


__all__ = [
    "FormalAxisPicardFamilyState",
    "formal_axis_picard_family_state",
]
