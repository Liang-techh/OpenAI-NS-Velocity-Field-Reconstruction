"""Fixed dyadic geometry from Eqs. (6.1)-(6.6).

Only the geometry is implemented. This does not select slow supports, separated
rectangles, transported waves, rounded frequencies or amplitudes. Q is a fixed
chart parameter; it must not be differentiated as the pointwise scale q(z,t).
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from .coordinates import validate_h

LAMBDA_G = 4 - math.sqrt(2)
T_G = 4 + math.sqrt(2)
B_G = math.sqrt(2) - 1
V_R = np.array([1.0, -B_G])
V_T = np.array([B_G, 1.0])
J_G = np.array([[3, 1], [1, 5]], dtype=object)  # Python integers avoid overflow
for _constant in (V_R, V_T, J_G):
    _constant.setflags(write=False)
RHO_G = math.log(LAMBDA_G) / math.log(T_G)
KAPPA_S = 1e-5


@dataclass(frozen=True)
class DyadicChart:
    ell: int
    h: float = 0.005

    def __post_init__(self):
        validate_h(self.h)
        if isinstance(self.ell, bool) or not isinstance(self.ell, int) or not 1 <= self.ell <= 1000:
            raise ValueError("binary64 chart implementation requires integer 1<=ell<=1000")

    @property
    def Q(self):
        return math.ldexp(1.0, -self.ell)

    @property
    def epsilon(self):
        return self.Q**self.h

    @property
    def S_star(self):
        return self.ell**2

    @property
    def d_r(self):
        return 2 * ((1 + self.h) * RHO_G - self.h * KAPPA_S)

    @property
    def covering_index(self):
        i = math.floor(((1 + self.h) * self.ell * math.log(2)
                        - 2 * math.log(self.ell)) / math.log(T_G))
        if i < 0:
            raise ValueError("paper covering requires a sufficiently large band with i>=0")
        return i

    @property
    def c_i(self):
        return math.exp(self.covering_index * math.log(T_G)
                        - (1 + self.h) * self.ell * math.log(2))

    @property
    def M_i(self):
        return math.exp(self.covering_index * math.log(LAMBDA_G)
                        - self.d_r / 2 * self.ell * math.log(2))

    def covering_matrix(self):
        """Exact integer J_g^i, NOT a guarantee of accurate floating fast phases."""
        return np.linalg.matrix_power(J_G, self.covering_index)

    def from_physical_tau(self, r: float, z: float, tau: float) -> tuple[float, float, float]:
        if (not all(math.isfinite(v) for v in (r, z, tau)) or r < 0 or tau < 0):
            raise ValueError("need finite r>=0, tau>=0, z")
        return r / math.sqrt(self.Q), z / self.Q**(0.5 - self.h), tau / self.Q

    def to_physical_tau(self, R: float, Z: float, T: float) -> tuple[float, float, float]:
        if not all(math.isfinite(v) for v in (R, Z, T)) or R < 0 or T < 0:
            raise ValueError("need finite R>=0, T>=0, Z")
        return R * math.sqrt(self.Q), Z * self.Q**(0.5 - self.h), T * self.Q

    def absolute_torus(self, r: float, t: float) -> np.ndarray:
        """(6.3), evaluated in binary64; not a high-frequency precision certificate."""
        if not math.isfinite(r) or r < 0 or not math.isfinite(t):
            raise ValueError("need finite r>=0 and t")
        return np.mod(V_R * r**self.d_r + V_T * t, 1.0)

    def normalize(self, value, kind: str):
        """Chart representatives: Q^A u, Q^(2A) p, Q^(2A+1/2) residual."""
        A = 0.5 + self.h
        exponents = {"velocity": A, "pressure": 2 * A, "residual": 2 * A + 0.5}
        if kind not in exponents:
            raise ValueError("kind must be velocity, pressure, or residual")
        factor = self.Q**exponents[kind]
        a = np.asarray(value, dtype=float)
        if factor == 0 or not np.all(np.isfinite(a)):
            raise ArithmeticError("normalization is outside binary64 range")
        return factor * a

    def derivative_coefficients(self, R: float) -> dict[str, float]:
        """Coefficients of t*, D_r, D_z, D_theta in (6.6), R>0."""
        if not math.isfinite(R) or R <= 0:
            raise ValueError("chart derivative coefficients require R>0")
        return {"time_slow": -self.epsilon, "time_fast": self.c_i,
                "radial_slow": 1.0, "radial_fast": self.M_i * self.d_r * R**(self.d_r - 1),
                "axial_slow": self.epsilon, "angular": 1 / R}
