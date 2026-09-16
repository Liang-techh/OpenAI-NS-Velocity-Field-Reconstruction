"""Clean-room exact checks for the public Kokuno stage-input covariance formulas.

This module intentionally reimplements only a small algebraic seam.  It does not
copy or execute the source repository's checker code and it does not certify a
complete stage construction.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

HALF = Fraction(1, 2)
ZERO = Fraction(0, 1)
ONE = Fraction(1, 1)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_PUBLIC_MAIN = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
SOURCE_RECORD_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_COMPONENT_PATH = "proof_sources/stage_inputs_audit/stage_inputs_body.tex"
SOURCE_COMPONENT_SHA256 = (
    "61aa4684591b964321b00304d58ca7e8e2bbf538cfc930bc8a6be219755bd7e3"
)
SOURCE_COVARIANCE_RESULT_SHA256 = (
    "31f2a1ae9bfff6cfad82d40a1eb866cf9d8247ab0eace15d8737ca643255c99d"
)
SOURCE_BUNDLE_SHA256 = (
    "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
)


@dataclass(frozen=True)
class StageInputCovarianceReceipt:
    """Exact receipt for one covariance-scaling and assembly microcheck."""

    h: Fraction
    A: Fraction
    q_exponent: Fraction
    Q_exponent: Fraction
    partition_square_sum: Fraction
    scaling_relation_ok: bool
    partition_square_ok: bool
    verified: bool


def _require_fraction(name: str, value: Fraction) -> None:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be fractions.Fraction, got {type(value).__name__}")


def _exact_square_sum(weights: Iterable[Fraction]) -> Fraction:
    total = ZERO
    for index, weight in enumerate(weights):
        _require_fraction(f"partition_weights[{index}]", weight)
        total += weight * weight
    return total


def verify_stage_input_covariance(
    *,
    h: Fraction,
    A: Fraction,
    partition_weights: Iterable[Fraction],
) -> StageInputCovarianceReceipt:
    """Check the public stage-input scaling and partition-square identities exactly.

    The published scaling seam is

        Q^(-2A+h) (Q/q)^(A+1/2)
          = Q^(-A+h+1/2) q^(-A-1/2).

    Under ``A = 1/2 + h`` the Q exponent must vanish.  For a partitioned
    physical assembly, the local covariance weights combine exactly when
    ``sum_beta eta_beta^2 = 1``.  No floating tolerance is used.
    """

    _require_fraction("h", h)
    _require_fraction("A", A)

    q_exponent = -A - HALF
    Q_exponent = -A + h + HALF
    partition_square_sum = _exact_square_sum(partition_weights)

    scaling_relation_ok = A == HALF + h and Q_exponent == ZERO
    partition_square_ok = partition_square_sum == ONE

    return StageInputCovarianceReceipt(
        h=h,
        A=A,
        q_exponent=q_exponent,
        Q_exponent=Q_exponent,
        partition_square_sum=partition_square_sum,
        scaling_relation_ok=scaling_relation_ok,
        partition_square_ok=partition_square_ok,
        verified=scaling_relation_ok and partition_square_ok,
    )
