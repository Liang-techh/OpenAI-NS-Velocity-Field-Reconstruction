"""Microbenchmark only: this does not measure full reconstruction throughput."""
import json
import math
import timeit
from statistics import median
import numpy as np
from openai_ns_reconstruction import LeadingProfile


def main():
    X, eta = 2.0, 0.7
    def U(x, e):
        return e * math.exp(-x)
    profile = LeadingProfile(E=lambda x, e: 0.0, U=U,
                             dU_deta=lambda x, e: math.exp(-x))
    def old_average():
        xs = np.linspace(0.0, X, 801)
        vals = np.array([U(float(x), eta) for x in xs])
        return float(np.sum((vals[1:] + vals[:-1]) * np.diff(xs) / 2) / X)
    def new_average():
        return profile.radial_average_U(X, eta)
    new_average()  # warm the cached rule
    old_seconds = median(timeit.repeat(old_average, number=1000, repeat=5))
    new_seconds = median(timeit.repeat(new_average, number=1000, repeat=5))
    exact = eta * -math.expm1(-X) / X
    print(json.dumps({"scope": "single smooth radial-average microbenchmark only",
                      "calls_per_repeat": 1000, "repeats": 5,
                      "old_seconds": old_seconds, "new_seconds": new_seconds,
                      "speedup": old_seconds / new_seconds,
                      "old_absolute_error": abs(old_average()-exact),
                      "new_absolute_error": abs(new_average()-exact)}, indent=2))


if __name__ == "__main__":
    main()
