"""Bounded pressure-only variable projection for constrained velocity candidates.

This diagnostic freezes the candidate velocity and either uses no forcing or a
caller-supplied *fixed* forcing vector. It then asks how much of the sampled
momentum residual can be removed by coefficients in one declared pressure-
gradient basis, subject to explicit coefficient bounds.

There is intentionally no optimizable forcing basis and no ``f = R`` escape
hatch. A low projected residual is a same-sample pressure-capacity result, not
independent PDE validation and not evidence that the velocity field matches any
particular public visualization.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import lsq_linear


@dataclass(frozen=True)
class PressureProjectionResult:
    parameter_labels: tuple[str, ...]
    coefficients: tuple[float, ...]
    active_lower_bounds: tuple[str, ...]
    active_upper_bounds: tuple[str, ...]
    design_rank: int
    design_nullity: int
    design_condition: float | None
    residual_rms_before: float
    residual_rms_after: float
    residual_max_before: float
    residual_max_after: float
    recoverable_fraction_rms: float
    forcing_mode: str
    solver_success: bool
    solver_status: int
    solver_message: str
    function_evaluations: int | None
    velocity_changed: bool = False
    independent_validation_required: bool = True
    truth_boundary: str = (
        "same-sample bounded pressure-capacity diagnostic only; "
        "not optimizer generalization, not independent PDE validation, "
        "not visual correspondence, and not identification of an OpenAI field"
    )

    def to_dict(self) -> dict:
        return asdict(self)


def _finite_vector(name: str, value, *, length: int | None = None) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a nonempty 1D vector")
    if length is not None and arr.size != length:
        raise ValueError(f"{name} must have length {length}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _finite_matrix(name: str, value) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 2 or min(arr.shape) == 0:
        raise ValueError(f"{name} must be a nonempty 2D matrix")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _labels(labels: Sequence[str] | None, count: int) -> tuple[str, ...]:
    if labels is None:
        return tuple(f"pressure_{i}" for i in range(count))
    result = tuple(str(item) for item in labels)
    if len(result) != count or len(set(result)) != count or any(not item for item in result):
        raise ValueError(
            "parameter_labels must be unique, nonempty, and match pressure columns"
        )
    return result


def _rank_and_condition(matrix: np.ndarray, rank_rtol: float) -> tuple[int, float | None]:
    singular = np.linalg.svd(matrix, compute_uv=False)
    if singular.size == 0 or singular[0] <= 0.0:
        return 0, None
    threshold = rank_rtol * singular[0]
    rank = int(np.count_nonzero(singular > threshold))
    if rank == 0:
        return 0, None
    return rank, float(singular[0] / singular[rank - 1])


def bounded_pressure_projection(
    momentum_residual,
    pressure_gradient_design,
    *,
    lower_bounds,
    upper_bounds,
    residual_scales,
    fixed_forcing=None,
    parameter_labels: Sequence[str] | None = None,
    rank_rtol: float = 1e-10,
) -> PressureProjectionResult:
    """Compute the best declared bounded pressure correction on fixed samples.

    The convention is

    ``R(c) = momentum_residual - fixed_forcing + G @ c``.

    ``momentum_residual`` therefore represents the frozen velocity/base-pressure
    side before this pressure correction and before the optional immutable force.
    The routine optimizes only ``c``. It accepts no forcing basis and cannot fit
    forcing coefficients.

    ``residual_scales`` are explicit positive row scales used only to pose the
    least-squares problem. Reported RMS/max values remain in the caller's
    original residual units.
    """
    residual = _finite_vector("momentum_residual", momentum_residual)
    design = _finite_matrix("pressure_gradient_design", pressure_gradient_design)
    if design.shape[0] != residual.size:
        raise ValueError("pressure_gradient_design row count must match residual length")

    parameter_count = design.shape[1]
    lower = _finite_vector("lower_bounds", lower_bounds, length=parameter_count)
    upper = _finite_vector("upper_bounds", upper_bounds, length=parameter_count)
    scales = _finite_vector("residual_scales", residual_scales, length=residual.size)

    if np.any(lower >= upper):
        raise ValueError("every pressure lower bound must be below its upper bound")
    if np.any(scales <= 0.0):
        raise ValueError("residual_scales must be strictly positive")
    if not np.isfinite(rank_rtol) or not 0.0 < rank_rtol < 1.0:
        raise ValueError("rank_rtol must lie in (0,1)")

    labels = _labels(parameter_labels, parameter_count)

    if fixed_forcing is None:
        forcing = np.zeros_like(residual)
        forcing_mode = "none"
    else:
        forcing = _finite_vector("fixed_forcing", fixed_forcing, length=residual.size)
        forcing_mode = "fixed"

    base = residual - forcing
    scaled_design = design / scales[:, np.newaxis]
    scaled_target = -base / scales

    result = lsq_linear(
        scaled_design,
        scaled_target,
        bounds=(lower, upper),
        method="trf",
        lsmr_tol="auto",
    )
    if not result.success or not np.all(np.isfinite(result.x)):
        raise RuntimeError(
            "bounded pressure projection failed: "
            f"status={result.status}, message={result.message}"
        )

    coeff = np.asarray(result.x, dtype=float)
    after = base + design @ coeff

    rank, condition = _rank_and_condition(scaled_design, rank_rtol)
    before_rms = float(np.sqrt(np.mean(base * base)))
    after_rms = float(np.sqrt(np.mean(after * after)))
    before_max = float(np.max(np.abs(base)))
    after_max = float(np.max(np.abs(after)))

    if before_rms <= np.finfo(float).eps:
        recovered = 0.0 if after_rms <= np.finfo(float).eps else float("-inf")
    else:
        recovered = float(1.0 - after_rms / before_rms)

    bound_scale = np.maximum(1.0, np.maximum(np.abs(lower), np.abs(upper)))
    bound_tol = 256.0 * np.finfo(float).eps * bound_scale
    active_lower = tuple(
        labels[i]
        for i in range(parameter_count)
        if abs(coeff[i] - lower[i]) <= bound_tol[i]
    )
    active_upper = tuple(
        labels[i]
        for i in range(parameter_count)
        if abs(coeff[i] - upper[i]) <= bound_tol[i]
    )

    nit = getattr(result, "nit", None)
    if nit is not None:
        nit = int(nit)

    return PressureProjectionResult(
        parameter_labels=labels,
        coefficients=tuple(float(value) for value in coeff),
        active_lower_bounds=active_lower,
        active_upper_bounds=active_upper,
        design_rank=rank,
        design_nullity=int(parameter_count - rank),
        design_condition=condition,
        residual_rms_before=before_rms,
        residual_rms_after=after_rms,
        residual_max_before=before_max,
        residual_max_after=after_max,
        recoverable_fraction_rms=recovered,
        forcing_mode=forcing_mode,
        solver_success=True,
        solver_status=int(result.status),
        solver_message=str(result.message),
        function_evaluations=nit,
    )
