"""Analytic Section 10 pressure-gradient bridge for one finite Section 9 prefix.

The pinned official spatial-localization source defines

    cutPressure(p)(t,x) = spatialCutoff(x) * p(t,x).

This module applies that exact algebra to an already-admitted
``Section9FinitePrefixJetCertificate`` at one caller-supplied spacetime point
and evaluates

    grad(c p) = c grad(p) + p grad(c)

with the repository's fixed Section 10 spatial cutoff and matching analytic
gradient.  Production code performs no finite differences.

This is deliberately only a bounded derivative adapter.  The finite Section 9
prefix remains provider-supplied, its base point is caller metadata, the
transition collar uses the repository's geometry-faithful executable cutoff
representative, and no actual Section 7/8 correction sequence, infinite
Eq. (9.21), endpoint extension, Navier--Stokes residual, forcing, energy
closure, blow-up closure, or paper-exact velocity is constructed here.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .section9_eq921_prefix_jet import Section9FinitePrefixJetCertificate
from .spatial_localization import (
    section10_spatial_cutoff,
    section10_spatial_cutoff_gradient,
)

_ZERO = (0, 0, 0, 0)
_DX = (0, 1, 0, 0)
_DY = (0, 0, 1, 0)
_DZ = (0, 0, 0, 1)


def _finite_scalar(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real scalar")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _finite_coordinate(value: object, name: str) -> float:
    return _finite_scalar(value, name)


def _readonly_vec3(value: object, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    frozen = np.array(out, dtype=float, copy=True)
    frozen.setflags(write=False)
    return frozen


@dataclass(frozen=True)
class Section9Section10LocalizedPressureGradientCertificate:
    """One-point analytic ``grad(c p)`` with explicit fail-closed truth flags."""

    point: tuple[float, float, float, float]
    q: object
    prefix_order: int
    derivative_order: int
    provider_id: str
    provider_revision: str
    provider_provenance: str
    pressure_value: float
    pressure_gradient: np.ndarray
    cutoff_value: float
    cutoff_gradient: np.ndarray
    localized_pressure: float
    localized_pressure_gradient: np.ndarray
    paper_section10_cutpressure_definition_applied: bool = True
    analytic_pressure_gradient_from_prefix_jet_verified: bool = True
    point_binding_to_actual_correction_field_verified: bool = False
    actual_section9_sequence_verified: bool = False
    eq_9_21_infinite_sum_constructed: bool = False
    section9_field_smooth_extension_through_t1_constructed: bool = False
    residual_artifact_ready: bool = False
    endpoint_residual_closure_verified: bool = False
    paper_exact_velocity_available: bool = False

    @property
    def formal_localized_pressure_gradient_ready(self) -> bool:
        return (
            self.derivative_order >= 1
            and self.paper_section10_cutpressure_definition_applied
            and self.analytic_pressure_gradient_from_prefix_jet_verified
            and not self.point_binding_to_actual_correction_field_verified
            and not self.actual_section9_sequence_verified
            and not self.eq_9_21_infinite_sum_constructed
            and not self.section9_field_smooth_extension_through_t1_constructed
            and not self.residual_artifact_ready
            and not self.endpoint_residual_closure_verified
            and not self.paper_exact_velocity_available
        )


def localize_section9_prefix_pressure_gradient(
    prefix: Section9FinitePrefixJetCertificate,
    *,
    x: float,
    y: float,
    z: float,
    t: float,
) -> Section9Section10LocalizedPressureGradientCertificate:
    """Evaluate the fixed Section 10 cut pressure and its spatial gradient.

    The cutoff is not caller-overridable.  ``prefix`` must pass its finite-prefix
    truth gate and contain first spatial derivatives of the weighted Section 9
    pressure prefix.
    """
    if not isinstance(prefix, Section9FinitePrefixJetCertificate):
        raise TypeError("prefix must be a Section9FinitePrefixJetCertificate")
    if not prefix.formal_prefix_jet_ready:
        raise ValueError("prefix must pass the finite-prefix jet truth gate")
    if prefix.derivative_order < 1:
        raise ValueError(
            "prefix derivative_order must be at least one to form grad(p)"
        )

    x0 = _finite_coordinate(x, "x")
    y0 = _finite_coordinate(y, "y")
    z0 = _finite_coordinate(z, "z")
    t0 = _finite_coordinate(t, "t")

    pressure_value = _finite_scalar(prefix.p_derivatives[_ZERO], "p")
    pressure_gradient = _readonly_vec3(
        np.array(
            [
                prefix.p_derivatives[_DX],
                prefix.p_derivatives[_DY],
                prefix.p_derivatives[_DZ],
            ],
            dtype=float,
        ),
        "grad(p)",
    )

    cutoff_value = _finite_scalar(
        section10_spatial_cutoff(x0, y0, z0, t0),
        "Section 10 cutoff",
    )
    cutoff_gradient = _readonly_vec3(
        section10_spatial_cutoff_gradient(x0, y0, z0, t0),
        "Section 10 cutoff gradient",
    )

    localized_pressure = _finite_scalar(
        cutoff_value * pressure_value,
        "Section 10 cut pressure",
    )
    localized_pressure_gradient = _readonly_vec3(
        cutoff_value * pressure_gradient + pressure_value * cutoff_gradient,
        "grad(Section 10 cut pressure)",
    )

    certificate = Section9Section10LocalizedPressureGradientCertificate(
        point=(x0, y0, z0, t0),
        q=prefix.q,
        prefix_order=prefix.prefix_order,
        derivative_order=prefix.derivative_order,
        provider_id=prefix.provider_id,
        provider_revision=prefix.provider_revision,
        provider_provenance=prefix.provider_provenance,
        pressure_value=pressure_value,
        pressure_gradient=pressure_gradient,
        cutoff_value=cutoff_value,
        cutoff_gradient=cutoff_gradient,
        localized_pressure=localized_pressure,
        localized_pressure_gradient=localized_pressure_gradient,
    )
    if not certificate.formal_localized_pressure_gradient_ready:
        raise ValueError("localized pressure-gradient certificate failed its truth gate")
    return certificate
