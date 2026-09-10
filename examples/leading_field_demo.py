"""Demonstrate the paper-exact kinematics with a toy smooth profile.

IMPORTANT: the profile below is NOT OpenAI's constructed profile.  It is only a smoke-test
for the coordinate map and velocity evaluator.  The repository does not label it paper-exact.
"""

import math
import numpy as np

from openai_ns_reconstruction import LeadingProfile
from openai_ns_reconstruction.velocity import blowup_probe
from openai_ns_reconstruction.verify import loglog_slope


def E(X: float, eta: float) -> float:
    return math.sqrt(2.0 * X) * math.exp(-X) * (1.0 + 0.1 * eta * eta)


def U(X: float, eta: float) -> float:
    return eta * math.exp(-X)


def dU_deta(X: float, eta: float) -> float:
    return math.exp(-X)


profile = LeadingProfile(E=E, U=U, dU_deta=dU_deta, name="toy-smoke-test", paper_exact=False)
h = 0.005
X_in = 0.5

taus = np.logspace(-2, -7, 10)
values = []
for tau in taus:
    v = blowup_probe(float(tau), X_in, profile, h=h)
    values.append(v.u_theta)
    print(f"tau={tau:.2e}  q={v.q:.2e}  X={v.X:.6f}  u_theta={v.u_theta:.6e}")

slope = loglog_slope(taus, np.asarray(values))
print(f"\nmeasured log-log exponent: {slope:.6f}")
print(f"predicted leading exponent: {-0.5-h:.6f}")
