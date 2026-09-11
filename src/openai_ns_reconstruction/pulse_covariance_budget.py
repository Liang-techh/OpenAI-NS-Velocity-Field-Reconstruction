"""Section 7 pulse-integral error budget behind Eq. (7.28).

This module isolates the analytic averaging step between the actual covariance
integral (7.27) and the finite-dimensional perturbation layer already provided
by :mod:`stress_cone_perturbation`.

For a signed slot ``sigma in {-1,+1}``, the paper uses

``s(v) = sigma * (u_*/2 + u_* v/L_s)``

and Lemma 7.4 gives the normalized tangential direction, in the frozen
``(N,K)`` frame, as

``(c_* sqrt(1+s(v)^2), -s(v)) + r(v)``

with a small ratio/frame defect ``r``.  At ``v=L_s/2`` the ideal direction is

``(-A_c, -sigma u_*)``,  ``A_c=-c_* sqrt(1+u_*^2)>0``.

If an upstream analytic argument certifies

* ``|r(v)| <= ratio_error_bound`` on the pulse support, and
* the normalized weighted first moment
  ``E_w[|v-L_s/2|/L_s] <= normalized_first_moment_bound``,

then the normalized covariance-column error obeys

``|e_sigma| <= ratio_error_bound
             + u_* sqrt(1+c_*^2) normalized_first_moment_bound``.

The proof is just the mean-value theorem for
``g(s)=(c_*sqrt(1+s^2),-s)`` followed by the triangle inequality under the
positive Eq. (7.27) weight ``w=psi^2 x^2``.  The paper's Lemma 7.4
``O(S_*^-1)`` ratio error and Gaussian first-moment ``O(S_*^-1/2)`` therefore
give Eq. (7.28)'s ``O(S_*^-1/2)`` column error.

This remains ``formal-structure``: the bounds supplied here must come from the
actual paper background/pulse construction.  Finite samples, quadrature data,
or arbitrary caller constants are not promoted to paper-exact certificates.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

from .stress_cone import normal_magnitude


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _positive(value: float, name: str) -> float:
    value = _finite(value, name)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")
    return value


def _nonnegative(value: float, name: str) -> float:
    value = _finite(value, name)
    if value < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _up_nonnegative(value: float, name: str) -> float:
    value = _nonnegative(value, name)
    if value == 0.0:
        return 0.0
    out = math.nextafter(value, math.inf)
    if not math.isfinite(out):
        raise OverflowError(f"{name} overflowed")
    return out


def _mul_up(x: float, y: float, name: str) -> float:
    x = _nonnegative(x, f"{name} left factor")
    y = _nonnegative(y, f"{name} right factor")
    if x == 0.0 or y == 0.0:
        return 0.0
    return _up_nonnegative(x * y, name)


def _add_up(x: float, y: float, name: str) -> float:
    x = _nonnegative(x, f"{name} left summand")
    y = _nonnegative(y, f"{name} right summand")
    if x == 0.0 and y == 0.0:
        return 0.0
    return _up_nonnegative(x + y, name)


def _vector2(values: Iterable[float], name: str) -> tuple[float, float]:
    values = tuple(float(x) for x in values)
    if len(values) != 2:
        raise ValueError(f"{name} must have length 2")
    if not all(math.isfinite(x) for x in values):
        raise ValueError(f"{name} must contain finite values")
    return values[0], values[1]


@dataclass(frozen=True)
class SignedPulseCovarianceErrorBudget:
    """Certified-input implication for one signed Eq. (7.28) pulse column.

    ``sigma`` is the paper's signed pulse label.  ``c_star`` and ``u_star`` are
    the reference constants from Proposition 7.5.  The two error inputs are
    *analytic upper bounds* that must be established upstream; this class does
    not estimate them from samples.
    """

    sigma: int
    c_star: float
    u_star: float
    ratio_error_bound: float
    normalized_first_moment_bound: float

    def __post_init__(self) -> None:
        if self.sigma not in (-1, 1):
            raise ValueError("sigma must be -1 or +1")
        c_star = _finite(self.c_star, "c_star")
        if c_star >= 0.0:
            raise ValueError("Proposition 7.5 signed convention requires c_star<0")
        _positive(self.u_star, "u_star")
        _nonnegative(self.ratio_error_bound, "ratio_error_bound")
        _nonnegative(
            self.normalized_first_moment_bound,
            "normalized_first_moment_bound",
        )

    @property
    def a_c(self) -> float:
        """Positive normal magnitude ``A_c=-c_*sqrt(1+u_*^2)``."""
        return normal_magnitude(self.c_star, self.u_star)

    @property
    def reference_direction(self) -> tuple[float, float]:
        """Eq. (7.28) reference direction in frozen ``(N,K)`` coordinates."""
        return -self.a_c, -float(self.sigma) * self.u_star

    def ideal_direction(self, s: float) -> tuple[float, float]:
        """Return ``g(s)=(c_*sqrt(1+s^2),-s)`` in ``(N,K)`` coordinates."""
        s = _finite(s, "s")
        normal = self.c_star * math.hypot(1.0, s)
        transverse = -s
        if not math.isfinite(normal):
            raise OverflowError("ideal normal component overflowed")
        return normal, transverse

    @property
    def ideal_direction_lipschitz_upper(self) -> float:
        """Global bound ``|g'(s)| <= sqrt(1+c_*^2)``.

        Indeed ``|g'(s)|^2 = 1 + c_*^2 s^2/(1+s^2)``.
        """
        return _up_nonnegative(math.hypot(1.0, self.c_star), "ideal direction Lipschitz bound")

    @property
    def concentration_error_upper(self) -> float:
        """Weighted center-drift contribution to the normalized column error."""
        first = _mul_up(
            self.u_star,
            self.ideal_direction_lipschitz_upper,
            "center-drift coefficient",
        )
        return _mul_up(
            first,
            self.normalized_first_moment_bound,
            "concentration error",
        )

    @property
    def normalized_error_bound(self) -> float:
        """The theorem-shaped upper bound for ``|e_sigma|`` in Eq. (7.28)."""
        return _add_up(
            self.ratio_error_bound,
            self.concentration_error_upper,
            "normalized pulse-column error",
        )

    def error_vector_for(self, normalized_column: Iterable[float]) -> tuple[float, float]:
        """Check one supplied normalized column against this analytic budget.

        This is a downstream consistency check only.  Supplying a column does
        not prove that it is the paper's Eq. (7.27) integral.
        """
        column = _vector2(normalized_column, "normalized_column")
        reference = self.reference_direction
        error = (column[0] - reference[0], column[1] - reference[1])
        norm = math.hypot(*error)
        norm_upper = 0.0 if norm == 0.0 else _up_nonnegative(norm, "column error norm")
        if norm_upper > self.normalized_error_bound:
            raise ValueError("normalized column error exceeds analytic Eq. (7.28) budget")
        return error


def uniform_two_sign_error_bound(
    minus: SignedPulseCovarianceErrorBudget,
    plus: SignedPulseCovarianceErrorBudget,
) -> float:
    """Return one ``delta`` suitable for the two-column perturbation layer.

    The two certificates must describe the same ``c_*`` and ``u_*`` and carry
    the expected signs.  Their maximum error bound can then be passed to
    ``CovariancePerturbationCertificate.normalized_error_bound`` once the actual
    Eq. (7.27) normalized columns have been constructed independently.
    """
    if minus.sigma != -1 or plus.sigma != 1:
        raise ValueError("expected minus.sigma=-1 and plus.sigma=+1")
    if minus.c_star != plus.c_star or minus.u_star != plus.u_star:
        raise ValueError("signed budgets must use identical c_star and u_star")
    return max(minus.normalized_error_bound, plus.normalized_error_bound)


def sqrt_scale_error_envelope(
    *,
    s_star: float,
    c_star: float,
    u_star: float,
    ratio_constant: float,
    first_moment_constant: float,
) -> float:
    """Collapse the paper's component rates to an explicit ``C/sqrt(S_*)``.

    If ``S_*>=1``, ``ratio_error <= ratio_constant/S_*`` and
    ``first_moment <= first_moment_constant/sqrt(S_*)``, then

    ``|e_sigma| <= (ratio_constant
                    + u_*sqrt(1+c_*^2) first_moment_constant)/sqrt(S_*)``.

    This function proves only the scalar implication.  It does not establish
    either component estimate for the real pulse.
    """
    s_star = _positive(s_star, "s_star")
    if s_star < 1.0:
        raise ValueError("sqrt-rate envelope requires s_star>=1")
    c_star = _finite(c_star, "c_star")
    if c_star >= 0.0:
        raise ValueError("Proposition 7.5 signed convention requires c_star<0")
    u_star = _positive(u_star, "u_star")
    ratio_constant = _nonnegative(ratio_constant, "ratio_constant")
    first_moment_constant = _nonnegative(first_moment_constant, "first_moment_constant")

    lipschitz = _up_nonnegative(math.hypot(1.0, c_star), "ideal direction Lipschitz bound")
    concentration_constant = _mul_up(
        _mul_up(u_star, lipschitz, "sqrt-rate drift coefficient"),
        first_moment_constant,
        "sqrt-rate concentration constant",
    )
    numerator = _add_up(ratio_constant, concentration_constant, "sqrt-rate numerator")
    return _up_nonnegative(numerator / math.sqrt(s_star), "sqrt-rate error envelope")
