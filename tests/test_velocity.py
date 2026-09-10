import math
import numpy as np

from openai_ns_reconstruction import LeadingProfile
from openai_ns_reconstruction.velocity import blowup_probe, leading_velocity_cartesian
from openai_ns_reconstruction.verify import loglog_slope


def toy_profile() -> LeadingProfile:
    def E(X, eta):
        return math.sqrt(2.0 * X) * math.exp(-X) * (1.0 + 0.1 * eta * eta)

    def U(X, eta):
        return eta * math.exp(-X)

    def dU_deta(X, eta):
        return math.exp(-X)

    return LeadingProfile(E=E, U=U, dU_deta=dU_deta, name="toy", paper_exact=False)


def test_cartesian_axis_is_finite():
    u = leading_velocity_cartesian(0.0, 0.0, 0.0, 0.8, toy_profile(), h=0.005)
    assert np.all(np.isfinite(u))
    assert u[0] == 0.0
    assert u[1] == 0.0


def test_probe_has_predicted_leading_exponent():
    h = 0.005
    taus = np.logspace(-2, -6, 8)
    vals = np.array([blowup_probe(float(t), 0.5, toy_profile(), h=h).u_theta for t in taus])
    slope = loglog_slope(taus, vals)
    assert abs(slope - (-0.5 - h)) < 1e-8
