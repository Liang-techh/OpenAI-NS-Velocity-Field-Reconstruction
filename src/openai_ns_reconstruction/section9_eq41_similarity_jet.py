"""Machine-derived finite spacetime jets for the Section 9 similarity coordinate.

This module removes one provider-input boundary from the Section 9 cutoff path.
Starting from the paper's Eq. (4.1),

    q - z^2 q^(2h) = tau,   tau = 1 - t,

it derives every derivative of ``q(t,x,y,z)`` through a requested finite total
order.  The implementation solves the implicit equation once with the existing
Eq. (4.1) coordinate solver and then uses a normalized bivariate Taylor
recurrence in ``(dt,dz)``.  No finite differences, interpolation, fitting, or
sampled derivative estimates are used.

This is still only an analytic derivative adapter.  It does not identify the
one-dimensional cutoff derivative table with the manuscript's fixed cutoff,
does not construct the infinite Eq. (9.21) sum or its t=1 extension, and does
not upgrade any velocity/forcing paper-exact claim.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral
from typing import Mapping

from .coordinates import solve_q_from_tau, validate_h
from .section9_cutoff_composition_jet import (
    Section9DerivedCutoffWeightJetCertificate,
    Section9ScalarCutoffDerivativeJet,
    Section9SimilarityCoordinateJet,
    derive_section9_cutoff_weight_jet,
)
from .section9_correction_extension_admission import (
    Section9CorrectionStageAdmissionCertificate,
)
from .section9_eq921_prefix_jet import MultiIndex4, required_spacetime_multiindices
from .section9_finite_prefix import Section9FinitePrefixEvaluationCertificate


_BiIndex = tuple[int, int]
_ZERO2: _BiIndex = (0, 0)


def _nonempty(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be nonempty")
    return value.strip()


def _order(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError("derivative_order must be a nonnegative integer")
    return int(value)


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real scalar")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _multiply(
    left: Mapping[_BiIndex, float],
    right: Mapping[_BiIndex, float],
    *,
    order: int,
) -> dict[_BiIndex, float]:
    out: dict[_BiIndex, float] = {}
    for (it, iz), left_value in left.items():
        for (jt, jz), right_value in right.items():
            key = (it + jt, iz + jz)
            if sum(key) <= order:
                value = out.get(key, 0.0) + left_value * right_value
                if not math.isfinite(value):
                    raise ArithmeticError("Eq. (4.1) Taylor multiplication overflowed")
                out[key] = value
    return out


def _fractional_power_series(
    q_series: Mapping[_BiIndex, float],
    *,
    q0: float,
    exponent: float,
    order: int,
) -> dict[_BiIndex, float]:
    """Return normalized Taylor coefficients of ``Q**exponent``.

    Because q0>0, write Q=q0*(1+R) and use the finite generalized-binomial
    expansion.  R has zero constant coefficient, so terms R^k with k>order
    cannot affect the requested total-order jet.
    """
    r = {
        key: value / q0
        for key, value in q_series.items()
        if key != _ZERO2 and value != 0.0
    }
    q0_power = q0**exponent
    if not math.isfinite(q0_power):
        raise ArithmeticError("Eq. (4.1) fractional power overflowed")

    result: dict[_BiIndex, float] = {_ZERO2: q0_power}
    power: dict[_BiIndex, float] = {_ZERO2: 1.0}
    binomial = 1.0
    for k in range(1, order + 1):
        power = _multiply(power, r, order=order)
        binomial *= (exponent - (k - 1)) / k
        factor = q0_power * binomial
        for key, value in power.items():
            updated = result.get(key, 0.0) + factor * value
            if not math.isfinite(updated):
                raise ArithmeticError("Eq. (4.1) fractional-power series overflowed")
            result[key] = updated
    return result


def _implicit_residual_series(
    q_series: Mapping[_BiIndex, float],
    *,
    q0: float,
    z0: float,
    tau0: float,
    exponent: float,
    order: int,
) -> dict[_BiIndex, float]:
    q_power = _fractional_power_series(
        q_series, q0=q0, exponent=exponent, order=order
    )
    z_squared: dict[_BiIndex, float] = {
        (0, 0): z0 * z0,
        (0, 1): 2.0 * z0,
        (0, 2): 1.0,
    }
    product = _multiply(z_squared, q_power, order=order)
    residual = dict(q_series)
    for key, value in product.items():
        residual[key] = residual.get(key, 0.0) - value

    # tau(t0+dt)=tau0-dt, hence -tau contributes -tau0 + dt.
    residual[_ZERO2] = residual.get(_ZERO2, 0.0) - tau0
    residual[(1, 0)] = residual.get((1, 0), 0.0) + 1.0
    return residual


@dataclass(frozen=True)
class Section9Eq41SimilarityCoordinateJetCertificate:
    """Finite q-jet derived directly from the Eq. (4.1) implicit relation."""

    q_jet: Section9SimilarityCoordinateJet
    tau: float
    z: float
    h: float
    implicit_jacobian: float
    relation_residual: float
    normalized_taylor_implicit_recurrence_verified: bool = True
    similarity_coordinate_jet_machine_derived_from_eq_4_1: bool = True
    used_finite_differences: bool = False
    paper_fixed_cutoff_derivatives_machine_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_eq41_jet_ready(self) -> bool:
        return (
            self.normalized_taylor_implicit_recurrence_verified
            and self.similarity_coordinate_jet_machine_derived_from_eq_4_1
            and self.implicit_jacobian > 0.0
            and math.isfinite(self.relation_residual)
            and not self.used_finite_differences
            and not self.paper_fixed_cutoff_derivatives_machine_verified
            and not self.paper_exact_velocity_available
        )


@dataclass(frozen=True)
class Section9Eq41CutoffCompositionCertificate:
    """Existing cutoff composition fed by a machine-derived Eq. (4.1) q-jet."""

    eq41_q_jet: Section9Eq41SimilarityCoordinateJetCertificate
    composition: Section9DerivedCutoffWeightJetCertificate
    similarity_coordinate_jet_machine_derived_from_eq_4_1: bool = True
    paper_fixed_cutoff_derivatives_machine_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_composition_ready(self) -> bool:
        return (
            self.eq41_q_jet.formal_eq41_jet_ready
            and self.composition.formal_composition_ready
            and self.similarity_coordinate_jet_machine_derived_from_eq_4_1
            and not self.paper_fixed_cutoff_derivatives_machine_verified
            and not self.paper_exact_velocity_available
        )


def derive_section9_similarity_coordinate_jet_from_eq41(
    *,
    tau: float,
    z: float,
    h: float,
    derivative_order: int,
    provider_id: str,
    provider_revision: str,
    provider_provenance: str,
) -> Section9Eq41SimilarityCoordinateJetCertificate:
    """Derive a finite ``(t,x,y,z)`` q-jet from Eq. (4.1).

    Let ``Q(dt,dz)`` be the normalized Taylor series for q.  For each total
    degree, the as-yet unknown coefficient enters the implicit residual with
    multiplier

        F_q = 1 - 2h z^2 q^(2h-1).

    The physical branch has positive F_q.  Setting each same-degree unknown to
    minus the already-known residual coefficient divided by F_q recursively
    determines the complete finite jet.  Since Eq. (4.1) is independent of x
    and y, all derivatives containing x or y are identically zero.
    """
    order = _order(derivative_order)
    tau0 = _finite(tau, "tau")
    if tau0 <= 0.0:
        raise ValueError("tau must be positive; the adapter approaches t=1 from below")
    z0 = _finite(z, "z")
    h0 = validate_h(_finite(h, "h"))
    source_id = _nonempty(provider_id, "provider_id")
    source_revision = _nonempty(provider_revision, "provider_revision")
    source_provenance = _nonempty(provider_provenance, "provider_provenance")

    q0 = solve_q_from_tau(z0, tau0, h0)
    exponent = 2.0 * h0
    implicit_jacobian = 1.0 - exponent * z0 * z0 * q0 ** (exponent - 1.0)
    if not math.isfinite(implicit_jacobian) or implicit_jacobian <= 0.0:
        raise ArithmeticError("Eq. (4.1) physical-branch implicit Jacobian must be positive")

    base_residual = q0 - z0 * z0 * q0**exponent - tau0
    residual_scale = max(1.0, abs(q0), abs(tau0), abs(z0 * z0 * q0**exponent))
    if not math.isfinite(base_residual) or abs(base_residual) > 5e-12 * residual_scale:
        raise ArithmeticError("Eq. (4.1) solver result does not satisfy the implicit relation")

    normalized_q: dict[_BiIndex, float] = {_ZERO2: q0}
    for total_degree in range(1, order + 1):
        for t_degree in range(total_degree + 1):
            z_degree = total_degree - t_degree
            key = (t_degree, z_degree)
            normalized_q[key] = 0.0
            residual = _implicit_residual_series(
                normalized_q,
                q0=q0,
                z0=z0,
                tau0=tau0,
                exponent=exponent,
                order=order,
            ).get(key, 0.0)
            coefficient = -residual / implicit_jacobian
            if not math.isfinite(coefficient):
                raise ArithmeticError("Eq. (4.1) implicit Taylor recurrence overflowed")
            normalized_q[key] = coefficient

    derivatives: dict[MultiIndex4, float] = {}
    for alpha in required_spacetime_multiindices(order):
        dt, dx, dy, dz = alpha
        if dx or dy:
            value = 0.0
        else:
            value = (
                normalized_q[(dt, dz)]
                * math.factorial(dt)
                * math.factorial(dz)
            )
        if not math.isfinite(value):
            raise ArithmeticError("derived Eq. (4.1) q derivative overflowed")
        derivatives[alpha] = float(value)

    q_jet = Section9SimilarityCoordinateJet(
        q=q0,
        derivative_order=order,
        derivatives=derivatives,
        provider_id=f"eq4.1-analytic:{source_id}",
        provider_revision=source_revision,
        provider_provenance=source_provenance,
    )
    certificate = Section9Eq41SimilarityCoordinateJetCertificate(
        q_jet=q_jet,
        tau=tau0,
        z=z0,
        h=h0,
        implicit_jacobian=implicit_jacobian,
        relation_residual=base_residual,
    )
    if not certificate.formal_eq41_jet_ready:
        raise ArithmeticError("Eq. (4.1) similarity-coordinate truth-status invariant failed")
    return certificate


def derive_section9_cutoff_weight_jet_from_eq41(
    evaluation: Section9FinitePrefixEvaluationCertificate,
    admission: Section9CorrectionStageAdmissionCertificate,
    cutoff_jet: Section9ScalarCutoffDerivativeJet,
    *,
    tau: float,
    z: float,
    h: float,
    provider_id: str,
    provider_revision: str,
    provider_provenance: str,
) -> Section9Eq41CutoffCompositionCertificate:
    """Feed a machine-derived Eq. (4.1) q-jet into the existing chi(a_j q) path."""
    eq41_q_jet = derive_section9_similarity_coordinate_jet_from_eq41(
        tau=tau,
        z=z,
        h=h,
        derivative_order=cutoff_jet.derivative_order,
        provider_id=provider_id,
        provider_revision=provider_revision,
        provider_provenance=provider_provenance,
    )
    composition = derive_section9_cutoff_weight_jet(
        evaluation, admission, eq41_q_jet.q_jet, cutoff_jet
    )
    certificate = Section9Eq41CutoffCompositionCertificate(
        eq41_q_jet=eq41_q_jet,
        composition=composition,
    )
    if not certificate.formal_composition_ready:
        raise ArithmeticError("Eq. (4.1) cutoff-composition truth-status invariant failed")
    return certificate
