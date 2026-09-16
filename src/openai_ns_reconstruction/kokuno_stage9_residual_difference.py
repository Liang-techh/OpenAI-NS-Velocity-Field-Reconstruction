"""Clean-room exact replay of the Stage-9 zero-order residual-difference identity.

This module is independently authored from the short public Cartesian formulas
in the Kokuno Navier--Stokes reader.  It intentionally does not copy the source
checker implementation.  The scope is the I=0 case of the public A38 expansion
and its scalar A39 majorant.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence, Tuple

Rat = Fraction
Vec3 = Tuple[Rat, Rat, Rat]
Mat3 = Tuple[Vec3, Vec3, Vec3]


def _rat(value: object, *, name: str) -> Rat:
    if isinstance(value, bool) or not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction")
    return value


def _vec3(values: Sequence[Rat], *, name: str) -> Vec3:
    if len(values) != 3:
        raise ValueError(f"{name} must have length 3")
    return tuple(_rat(v, name=f"{name}[{i}]") for i, v in enumerate(values))  # type: ignore[return-value]


def _mat3(values: Sequence[Sequence[Rat]], *, name: str) -> Mat3:
    if len(values) != 3:
        raise ValueError(f"{name} must have 3 rows")
    return tuple(_vec3(row, name=f"{name}[{i}]") for i, row in enumerate(values))  # type: ignore[return-value]


def _add_vec(a: Vec3, b: Vec3) -> Vec3:
    return tuple(a[i] + b[i] for i in range(3))  # type: ignore[return-value]


def _add_mat(a: Mat3, b: Mat3) -> Mat3:
    return tuple(_add_vec(a[i], b[i]) for i in range(3))  # type: ignore[return-value]


def _sub_vec(a: Vec3, b: Vec3) -> Vec3:
    return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]


def ns_residual_zero_order(
    *,
    velocity: Sequence[Rat],
    dt_velocity: Sequence[Rat],
    grad_velocity: Sequence[Sequence[Rat]],
    lap_velocity: Sequence[Rat],
    grad_pressure: Sequence[Rat],
) -> Vec3:
    """Return the viscosity-one Cartesian Navier--Stokes residual at one jet.

    ``grad_velocity[i][k]`` denotes ``partial_k velocity_i`` and
    ``lap_velocity[i]`` denotes ``sum_k partial_k^2 velocity_i``.
    """

    u = _vec3(velocity, name="velocity")
    dt_u = _vec3(dt_velocity, name="dt_velocity")
    grad_u = _mat3(grad_velocity, name="grad_velocity")
    lap_u = _vec3(lap_velocity, name="lap_velocity")
    grad_p = _vec3(grad_pressure, name="grad_pressure")

    out = []
    for i in range(3):
        transport = sum(u[k] * grad_u[i][k] for k in range(3))
        out.append(dt_u[i] + transport - lap_u[i] + grad_p[i])
    return tuple(out)  # type: ignore[return-value]


def expanded_difference_zero_order(
    *,
    v: Sequence[Rat],
    e: Sequence[Rat],
    dt_e: Sequence[Rat],
    grad_v: Sequence[Sequence[Rat]],
    grad_e: Sequence[Sequence[Rat]],
    lap_e: Sequence[Rat],
    grad_pi: Sequence[Rat],
) -> Vec3:
    """Return the complete I=0 A38 residual-difference expansion.

    Every linear, cross, and quadratic term is retained:
    ``e_t - Delta e + grad pi + v.grad(e) + e.grad(v) + e.grad(e)``.
    """

    vv = _vec3(v, name="v")
    ee = _vec3(e, name="e")
    dte = _vec3(dt_e, name="dt_e")
    gv = _mat3(grad_v, name="grad_v")
    ge = _mat3(grad_e, name="grad_e")
    le = _vec3(lap_e, name="lap_e")
    gpi = _vec3(grad_pi, name="grad_pi")

    out = []
    for i in range(3):
        interaction = sum(
            vv[k] * ge[i][k] + ee[k] * gv[i][k] + ee[k] * ge[i][k]
            for k in range(3)
        )
        out.append(dte[i] - le[i] + gpi[i] + interaction)
    return tuple(out)  # type: ignore[return-value]


@dataclass(frozen=True)
class ResidualDifferenceReplay:
    direct: Vec3
    expanded: Vec3
    residual: Vec3

    @property
    def exact(self) -> bool:
        return self.residual == (Rat(0), Rat(0), Rat(0))


def replay_zero_order_difference(
    *,
    v: Sequence[Rat],
    e: Sequence[Rat],
    dt_v: Sequence[Rat],
    dt_e: Sequence[Rat],
    grad_v: Sequence[Sequence[Rat]],
    grad_e: Sequence[Sequence[Rat]],
    lap_v: Sequence[Rat],
    lap_e: Sequence[Rat],
    grad_p: Sequence[Rat],
    grad_pi: Sequence[Rat],
) -> ResidualDifferenceReplay:
    """Compare the direct nonlinear residual difference with the A38 expansion."""

    vv = _vec3(v, name="v")
    ee = _vec3(e, name="e")
    dtv = _vec3(dt_v, name="dt_v")
    dte = _vec3(dt_e, name="dt_e")
    gv = _mat3(grad_v, name="grad_v")
    ge = _mat3(grad_e, name="grad_e")
    lv = _vec3(lap_v, name="lap_v")
    le = _vec3(lap_e, name="lap_e")
    gp = _vec3(grad_p, name="grad_p")
    gpi = _vec3(grad_pi, name="grad_pi")

    direct_full = ns_residual_zero_order(
        velocity=_add_vec(vv, ee),
        dt_velocity=_add_vec(dtv, dte),
        grad_velocity=_add_mat(gv, ge),
        lap_velocity=_add_vec(lv, le),
        grad_pressure=_add_vec(gp, gpi),
    )
    base = ns_residual_zero_order(
        velocity=vv,
        dt_velocity=dtv,
        grad_velocity=gv,
        lap_velocity=lv,
        grad_pressure=gp,
    )
    direct = _sub_vec(direct_full, base)
    expanded = expanded_difference_zero_order(
        v=vv,
        e=ee,
        dt_e=dte,
        grad_v=gv,
        grad_e=ge,
        lap_e=le,
        grad_pi=gpi,
    )
    return ResidualDifferenceReplay(
        direct=direct,
        expanded=expanded,
        residual=_sub_vec(direct, expanded),
    )


def componentwise_a39_majorant_m0(*, v2_bound: Rat, e2_bound: Rat) -> Rat:
    """Return the exact scalar m=0 majorant from the A39 counting argument.

    The public expansion has five linear derivative/pressure slots, six cross
    transport products, and three quadratic error products.  With supplied
    certified nonnegative bounds ``|V|_2 <= v2_bound`` and
    ``|e|_2 <= e2_bound``, each scalar residual component is bounded by

        (5 + 6*v2_bound + 3*e2_bound) * e2_bound.

    This is the rational scalar bound before the optional Euclidean-output
    ``sqrt(3)`` factor, so no floating approximation to sqrt(3) is introduced.
    """

    vb = _rat(v2_bound, name="v2_bound")
    eb = _rat(e2_bound, name="e2_bound")
    if vb < 0 or eb < 0:
        raise ValueError("jet bounds must be nonnegative")
    return (Rat(5) + Rat(6) * vb + Rat(3) * eb) * eb


def verify_componentwise_majorant(*, difference: Sequence[Rat], majorant: Rat) -> bool:
    diff = _vec3(difference, name="difference")
    bound = _rat(majorant, name="majorant")
    if bound < 0:
        raise ValueError("majorant must be nonnegative")
    return all(abs(x) <= bound for x in diff)
