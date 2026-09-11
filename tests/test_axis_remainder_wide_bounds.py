from decimal import Decimal
import math
import sys

import pytest

from openai_ns_reconstruction.axis_componentwise_input_bounds import (
    componentwise_analytic_input_norm_certificate,
)
from openai_ns_reconstruction.axis_remainder_bounds import (
    AxisDataNormBounds,
    NaturalOperatorNormBounds,
    natural_remainder_bound_certificate,
)
from openai_ns_reconstruction.axis_remainder_wide_bounds import (
    natural_remainder_bound_certificate_wide,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.schedule_analytic_neighborhood import (
    certify_actual_schedule_analytic_inputs,
)


def _axis_data() -> AxisDataNormBounds:
    return AxisDataNormBounds(
        A=0.3,
        D=0.2,
        h=0.01,
        one=1.0,
        eta=0.5,
        d=0.25,
        inverse_l=0.75,
        u_star=0.4,
        u_star_eta=0.6,
        w_star=0.7,
        h_star=0.8,
        normalized_gradient=0.9,
        z_star=1.1,
    )


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_wide_ledger_has_the_same_exact_simple_controlled_case() -> None:
    # Same deliberately simple Controlled bookkeeping oracle used by the
    # binary64 module: only the two dot terms survive, so the exact answer is 2.
    ops = NaturalOperatorNormBounds(
        product=1.0,
        average=0.0,
        primitive=0.0,
        parameter_primitive=0.0,
        mul_y=0.0,
        j1=0.0,
        j2=0.0,
        param1=0.0,
        param2=0.0,
        dot1=1.0,
        dot2=1.0,
        mixed1=0.0,
        mixed2=0.0,
        resolvent=1.0,
    )
    data = AxisDataNormBounds(
        A=0.0,
        D=0.0,
        h=0.0,
        one=0.0,
        eta=0.0,
        d=0.0,
        inverse_l=1.0,
        u_star=0.0,
        u_star_eta=0.0,
        w_star=1.0,
        h_star=0.0,
        normalized_gradient=0.0,
        z_star=0.0,
    )

    cert = natural_remainder_bound_certificate_wide(
        operators=ops,
        data=data,
        amplitude_norm_upper=0.0,
    )

    assert cert.reference_norm_upper == Decimal(0)
    assert cert.radius_upper == Decimal(1)
    assert cert.remainder_bound_upper == Decimal(2)
    assert cert.remainder_lipschitz_upper == Decimal(2)


def test_wide_ledger_cross_checks_binary64_when_binary64_is_representable() -> None:
    ops = NaturalOperatorNormBounds.pinned_axis_operators(
        epsilon=0.5,
        resolvent_norm_upper=7.0,
    )
    data = _axis_data()

    narrow = natural_remainder_bound_certificate(
        operators=ops,
        data=data,
        amplitude_norm_upper=0.2,
    )
    wide = natural_remainder_bound_certificate_wide(
        operators=ops,
        data=data,
        amplitude_norm_upper=0.2,
    )

    assert float(wide.reference_norm_upper) == pytest.approx(narrow.reference_norm_upper)
    assert float(wide.radius_upper) == pytest.approx(narrow.radius_upper)
    assert float(wide.remainder_bound_upper) == pytest.approx(
        narrow.remainder_bound_upper,
        rel=1e-12,
    )
    assert float(wide.remainder_lipschitz_upper) == pytest.approx(
        narrow.remainder_lipschitz_upper,
        rel=1e-12,
    )


def test_actual_schedule_remainder_majorant_survives_past_binary64_overflow() -> None:
    # This uses the already-landed actual SchedulePressure analytic certificate,
    # not hand-picked coefficient fields.  The old float propagation is known
    # to stop at this schedule because an intermediate remainder majorant is
    # non-finite.  Wide Decimal arithmetic must carry the same positive ledger
    # to a finite result rather than replacing it by a smaller surrogate.
    data = _schedule_data()
    schedule = certify_actual_schedule_analytic_inputs(data, 0.05)
    analytic = componentwise_analytic_input_norm_certificate(
        neighborhood_radius=schedule.neighborhood.radius,
        field_value_sup_upper=schedule.neighborhood.field_bounds,
    )
    operators = analytic.operator_norm_bounds()
    axis_data = analytic.axis_data_norm_bounds(h=data.h)

    wide = natural_remainder_bound_certificate_wide(
        operators=operators,
        data=axis_data,
        amplitude_norm_upper=analytic.amplitude_norm_upper,
    )

    assert wide.reference_norm_upper.is_finite()
    assert wide.radius_upper.is_finite()
    assert wide.remainder_bound_upper.is_finite()
    assert wide.remainder_lipschitz_upper.is_finite()
    assert wide.remainder_bound_upper > 0
    assert wide.remainder_lipschitz_upper > 0

    binary64_max = Decimal.from_float(sys.float_info.max)
    assert (
        wide.remainder_bound_upper > binary64_max
        or wide.remainder_lipschitz_upper > binary64_max
    )
    assert wide.decimal_order_bound is not None
    assert wide.decimal_order_lipschitz is not None


def test_wide_ledger_rejects_nonfinite_amplitude_input() -> None:
    with pytest.raises(ValueError, match="amplitude_norm_upper"):
        natural_remainder_bound_certificate_wide(
            operators=NaturalOperatorNormBounds.pinned_axis_operators(
                epsilon=1.0,
                resolvent_norm_upper=1.0,
            ),
            data=_axis_data(),
            amplitude_norm_upper=math.inf,
        )
