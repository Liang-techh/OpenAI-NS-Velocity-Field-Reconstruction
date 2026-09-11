from decimal import Decimal, localcontext
import math
import sys

import pytest

from openai_ns_reconstruction.natural_scale_selection import select_natural_scale
from openai_ns_reconstruction.natural_scale_selection_wide import (
    contraction_threshold_wide,
    select_natural_scale_wide,
    stability_scale_wide,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.stage1_scale_chain_wide import (
    diagnose_actual_schedule_scale_chain_wide,
)


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_wide_threshold_algebra_keeps_exact_positive_structure() -> None:
    bound = Decimal(2)
    lip = Decimal(3)

    assert contraction_threshold_wide(bound, lip) == Decimal(6)
    stability = stability_scale_wide(bound)

    # Independent rational oracle for 1 + (14000/9) * 2.  The production
    # value is rounded upward at finite precision, so it must enclose the same
    # high-precision rational evaluation rather than undershoot it.
    with localcontext() as ctx:
        ctx.prec = 90
        oracle = Decimal(1) + Decimal(28000) / Decimal(9)
    assert stability >= oracle
    assert stability - oracle < Decimal("1e-80")


def test_wide_selection_cross_checks_narrow_when_float_is_representable() -> None:
    narrow = select_natural_scale(
        remainder_bound=0.01,
        remainder_lipschitz=0.25,
        phase_real_part_sup=0.2,
    )
    wide = select_natural_scale_wide(
        remainder_bound=0.01,
        remainder_lipschitz=0.25,
        phase_real_part_sup=0.2,
    )

    assert float(wide.Lambda) == pytest.approx(narrow.Lambda, rel=1e-15)
    assert float(wide.C.exponent_upper) == pytest.approx(math.log(narrow.C), rel=1e-15)
    assert wide.admits_symbolic(wide.Lambda, wide.C.exponent_upper)
    assert not wide.admits_symbolic(
        wide.Lambda,
        wide.C.exponent_upper - Decimal("1e-20"),
    )


def test_symbolic_C_survives_beyond_any_direct_float_exponential() -> None:
    wide = select_natural_scale_wide(
        remainder_bound=Decimal("1e500"),
        remainder_lipschitz=Decimal("1e499"),
        phase_real_part_sup=Decimal("0.25"),
    )

    assert wide.Lambda.is_finite()
    assert wide.Lambda > Decimal.from_float(sys.float_info.max)
    assert wide.C.exponent_upper.is_finite()
    with localcontext() as ctx:
        ctx.prec = 90
        oracle_exponent = wide.Lambda * Decimal("0.25")
    assert wide.C.exponent_upper >= oracle_exponent
    assert (
        (wide.C.exponent_upper - oracle_exponent) / wide.C.exponent_upper
        < Decimal("1e-80")
    )
    # C itself is intentionally represented as exp(exponent_upper), so no
    # binary64/Decimal exponential overflow is needed to carry the theorem choice.
    assert wide.C.admits_exponential_exponent(wide.C.exponent_upper)


def test_actual_schedule_reaches_wide_lambda_C_without_surrogate_constants() -> None:
    diagnostic = diagnose_actual_schedule_scale_chain_wide(_schedule_data(), 0.05)

    assert diagnostic.upstream_obstruction is None
    assert diagnostic.reached_wide_remainder
    assert diagnostic.reached_wide_scale_selection
    assert diagnostic.remainder is not None
    assert diagnostic.scale is not None

    binary64_max = Decimal.from_float(sys.float_info.max)
    assert (
        diagnostic.remainder.remainder_bound_upper > binary64_max
        or diagnostic.remainder.remainder_lipschitz_upper > binary64_max
    )
    assert diagnostic.scale.Lambda.is_finite()
    assert diagnostic.scale.C.exponent_upper.is_finite()
    assert diagnostic.scale.C.exponent_upper >= 0
    assert diagnostic.scale.admits_symbolic(
        diagnostic.scale.Lambda,
        diagnostic.scale.C.exponent_upper,
    )


def test_wide_selection_rejects_uncertified_domains() -> None:
    with pytest.raises(ValueError, match="remainder_bound"):
        select_natural_scale_wide(
            remainder_bound=Decimal("-1"),
            remainder_lipschitz=Decimal(0),
            phase_real_part_sup=Decimal(0),
        )
    with pytest.raises(ValueError, match="phase_real_part_sup"):
        select_natural_scale_wide(
            remainder_bound=Decimal(0),
            remainder_lipschitz=Decimal(0),
            phase_real_part_sup=float("nan"),
        )
