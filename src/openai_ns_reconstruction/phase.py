"""Section 7.1 phase and tangent-frame formulas from the OpenAI NS paper.

This module implements the local algebra in Eqs. (7.1)--(7.8) for a single
slow box after the background jet and representative have been supplied by
the caller. It does *not* choose the paper's slow boxes, prove the uniform
estimates (7.9)--(7.11), solve the amplitude ODE, or claim a paper-exact wave.

Cross-check source: openai/NavierStokesAndEuler@f9e8bc5,
NavierStokes/PhaseCalculus.lean (phase, phaseNormal, backwardMaterialOp_phase).
The sign in :func:`backward_material_phase_defect` follows the paper/official
Lean convention ``t_* = partial_v - epsilon partial_T``.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class TangentialBaseJet:
    """Value and first slow derivatives of chart angular/axial speeds ``F,G``.

    Derivatives are with respect to the normalized slow coordinates ``R,Z,T``.
    """

    F: float
    G: float
    F_R: float = 0.0
    F_Z: float = 0.0
    F_T: float = 0.0
    G_R: float = 0.0
    G_Z: float = 0.0
    G_T: float = 0.0

    def __post_init__(self) -> None:
        if not all(math.isfinite(float(x)) for x in (
            self.F, self.G, self.F_R, self.F_Z, self.F_T,
            self.G_R, self.G_Z, self.G_T,
        )):
            raise ValueError("base jet entries must be finite")


@dataclass(frozen=True)
class ReferenceFrame:
    """Frozen tangential frame/growth data preceding Eq. (7.1)."""

    N: np.ndarray
    K: np.ndarray
    lambda0: float
    c0: float


@dataclass(frozen=True)
class PrimaryPhaseParameters:
    """Local phase choices from Eqs. (7.2)--(7.3)."""

    sigma: int
    epsilon: float
    k: int
    B_s: float
    u_star: float
    L_s: float
    p: float
    p_z: float
    x0: float
    angular_integer: int

    def s(self, v: float) -> float:
        if not math.isfinite(v):
            raise ValueError("v must be finite")
        return self.sigma * (self.u_star / 2.0 + self.u_star * v / self.L_s)


def _vec2(value, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (2,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite two-vector")
    return out


def _vec3(value, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be a finite three-vector")
    return out


def reference_frame(F0: float, g0) -> ReferenceFrame:
    """Construct ``N,K,lambda0,c0`` from the frozen base shear before (7.1)."""
    if not math.isfinite(F0):
        raise ValueError("F0 must be finite")
    g = _vec2(g0, "g0")
    g_norm = float(np.linalg.norm(g))
    if g_norm == 0.0:
        raise ValueError("g0 must be nonzero")
    N = g / g_norm
    # Tangential component order is (theta,z): K=N^perp=(-N_z,N_theta).
    K = np.array([-N[1], N[0]], dtype=float)
    lambda_sq = -2.0 * F0 * N[0] * (2.0 * F0 * N[0] + g_norm)
    denom = 2.0 * F0 * N[0]
    if not (lambda_sq > 0.0) or denom == 0.0:
        raise ValueError("profile data do not satisfy the positive-growth frame condition")
    lambda0 = math.sqrt(lambda_sq)
    c0 = lambda0 / denom
    if not c0 < 0.0:
        raise ValueError("paper frame requires c0<0")
    N.setflags(write=False)
    K.setflags(write=False)
    return ReferenceFrame(N=N, K=K, lambda0=lambda0, c0=c0)


def stress_cone_margins(target, frame: ReferenceFrame) -> tuple[float, float]:
    """Return the two strict margins corresponding to Eq. (7.1).

    Positive outputs mean ``T.N<0`` and ``|c0 (T.K)/(T.N)|<1``. The zero
    target is rejected because the paper interprets annulus edges by a one-sided
    limiting stress direction, not by division at the zero vector.
    """
    T = _vec2(target, "target")
    dot_n = float(T @ frame.N)
    if dot_n == 0.0:
        raise ValueError("Eq. (7.1) requires a nonzero N component / limiting edge direction")
    first = -dot_n
    second = 1.0 - abs(frame.c0 * float(T @ frame.K) / dot_n)
    return first, second


def _nearest_nonzero_integer(x: float) -> int:
    """Nearest nonzero integer with deterministic half-away-from-zero ties."""
    if not math.isfinite(x):
        raise ValueError("rounding target must be finite")
    if x >= 0.0:
        n = math.floor(x + 0.5)
    else:
        n = math.ceil(x - 0.5)
    if n == 0:
        # Among nonzero integers the nearest is sign(x); at x=0 the tie rule is +1.
        n = 1 if x >= 0.0 else -1
    return int(n)


def primary_phase_parameters(
    *, epsilon: float, lambda0: float, u_star: float, L_s: float,
    sigma: int, R0: float, K, g0,
) -> PrimaryPhaseParameters:
    """Instantiate the local frequency/phase parameters in Eq. (7.2).

    ``K`` and ``g0`` use tangential component order ``(theta,z)``. The paper
    leaves a fixed tie rule for rounding ``k*tilde(p)`` unspecified; this
    implementation uses half-away-from-zero ties and records the resulting
    exact nonzero angular integer.
    """
    if sigma not in (-1, 1):
        raise ValueError("sigma must be +1 or -1")
    if not all(math.isfinite(x) and x > 0.0 for x in (epsilon, lambda0, u_star, L_s, R0)):
        raise ValueError("epsilon, lambda0, u_star, L_s and R0 must be finite and positive")
    K2 = _vec2(K, "K")
    g = _vec2(g0, "g0")
    g2 = float(g @ g)
    if g2 == 0.0:
        raise ValueError("g0 must be nonzero")
    k = int(math.ceil(epsilon ** -0.5))
    B_s_sq = lambda0 / (epsilon * k * k * (1.0 + u_star * u_star) ** 1.5)
    B_s = math.sqrt(B_s_sq)
    unrounded = B_s * (K2 - sigma * u_star * g / (L_s * g2))
    p_tilde = R0 * float(unrounded[0])
    p_z = float(unrounded[1])
    angular_integer = _nearest_nonzero_integer(k * p_tilde)
    p = angular_integer / k
    x0 = sigma * B_s * u_star / 2.0
    return PrimaryPhaseParameters(
        sigma=sigma, epsilon=epsilon, k=k, B_s=B_s, u_star=u_star,
        L_s=L_s, p=p, p_z=p_z, x0=x0,
        angular_integer=angular_integer,
    )


def phase_value(
    R: float, Z: float, theta: float, v: float,
    params: PrimaryPhaseParameters, jet: TangentialBaseJet,
) -> float:
    """Equation (7.3): ``Phi=p theta + pz Z/eps + x0 R-v(pF+pzG)``."""
    if not all(math.isfinite(x) for x in (R, Z, theta, v)):
        raise ValueError("slot coordinates must be finite")
    return (
        params.p * theta + params.p_z * Z / params.epsilon + params.x0 * R
        - v * (params.p * jet.F + params.p_z * jet.G)
    )


def phase_normal(
    R: float, v: float, params: PrimaryPhaseParameters, jet: TangentialBaseJet,
) -> np.ndarray:
    """Equation (7.4): normalized cylindrical gradient ``n_Phi=grad_* Phi``."""
    if not math.isfinite(R) or R <= 0.0 or not math.isfinite(v):
        raise ValueError("Eq. (7.4) requires finite R>0 and v")
    p, pz, eps = params.p, params.p_z, params.epsilon
    return np.array([
        params.x0 - v * (p * jet.F_R + pz * jet.G_R),
        p / R,
        pz - eps * v * (p * jet.F_Z + pz * jet.G_Z),
    ], dtype=float)


def phase_normal_v_derivative(params: PrimaryPhaseParameters, jet: TangentialBaseJet) -> np.ndarray:
    """Exact ``partial_v n_Phi`` obtained by differentiating Eq. (7.4)."""
    p, pz, eps = params.p, params.p_z, params.epsilon
    return np.array([
        -(p * jet.F_R + pz * jet.G_R),
        0.0,
        -eps * (p * jet.F_Z + pz * jet.G_Z),
    ], dtype=float)


def backward_material_phase_defect(
    *, b: float, v: float, params: PrimaryPhaseParameters, jet: TangentialBaseJet,
) -> float:
    """Exact phase-material defect for the paper's backward slow-time operator.

    The operator is ``partial_v-eps*partial_T+b*partial_R+F*partial_theta
    +eps*G*partial_Z``. This matches the identity formalized in
    ``NavierStokes/PhaseCalculus.lean`` and is the algebra behind the
    ``E_ik`` estimate in Eq. (7.9).
    """
    if not math.isfinite(b) or not math.isfinite(v):
        raise ValueError("b and v must be finite")
    p, pz, eps = params.p, params.p_z, params.epsilon
    H_R = p * jet.F_R + pz * jet.G_R
    H_T = p * jet.F_T + pz * jet.G_T
    H_Z = p * jet.F_Z + pz * jet.G_Z
    return b * params.x0 - v * (b * H_R - eps * H_T + eps * jet.G * H_Z)


def shear_matrix(R: float, jet: TangentialBaseJet) -> np.ndarray:
    """The matrix ``K`` in Eq. (7.6)."""
    if not math.isfinite(R) or R <= 0.0:
        raise ValueError("Eq. (7.6) requires finite R>0")
    return np.array([
        [0.0, -2.0 * jet.F, 0.0],
        [2.0 * jet.F + R * jet.F_R, 0.0, 0.0],
        [jet.G_R, 0.0, 0.0],
    ], dtype=float)


def projected_evolution_operator(normal, normal_v_derivative, K) -> np.ndarray:
    """Projected amplitude operator ``A_Phi`` in Eq. (7.6)."""
    n = _vec3(normal, "normal")
    nv = _vec3(normal_v_derivative, "normal_v_derivative")
    mat = np.asarray(K, dtype=float)
    if mat.shape != (3, 3) or not np.all(np.isfinite(mat)):
        raise ValueError("K must be a finite 3x3 matrix")
    norm2 = float(n @ n)
    if norm2 == 0.0:
        raise ValueError("phase normal must be nonzero")
    row = n @ mat - nv
    return -mat + np.outer(n, row) / norm2


def tangent_frame(normal, c0: float) -> np.ndarray:
    """Return the real 3x2 frame ``B`` from Eqs. (7.7)--(7.8)."""
    n = _vec3(normal, "normal")
    if not math.isfinite(c0) or c0 == 0.0:
        raise ValueError("c0 must be finite and nonzero")
    tan = n[1:]
    tan_norm = float(np.linalg.norm(tan))
    if tan_norm == 0.0:
        raise ValueError("Eq. (7.7) requires nonzero tangential phase normal")
    K_a = tan / tan_norm
    N_a = np.array([K_a[1], -K_a[0]], dtype=float)
    s_a = n[0] / tan_norm
    U = np.array([
        [1.0, 0.0],
        [-s_a * K_a[0], N_a[0]],
        [-s_a * K_a[1], N_a[1]],
    ], dtype=float)
    scale = c0 * math.sqrt(1.0 + s_a * s_a)
    mixing = np.array([[1.0, 1.0], [scale, -scale]], dtype=float)
    return U @ mixing


def harmonic_is_single_valued(
    params: PrimaryPhaseParameters, harmonic: int = 1, *, atol: float = 1e-12,
) -> bool:
    """Check the angular-periodicity consequence of ``k*p in Z\\{0}``.

    This checks the carrier phase shift rather than evaluating a potentially
    ill-conditioned complex exponential at an arbitrary large base phase.
    """
    if isinstance(harmonic, bool) or not isinstance(harmonic, int):
        raise ValueError("harmonic must be an integer")
    shift_turns = params.k * harmonic * params.p
    return params.angular_integer != 0 and abs(shift_turns - round(shift_turns)) <= atol
