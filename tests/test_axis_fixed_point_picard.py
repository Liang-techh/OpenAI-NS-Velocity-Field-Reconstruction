from decimal import Decimal, localcontext

import pytest

from openai_ns_reconstruction.axis_fixed_point_picard import (
    NaturalPicardContractionCertificate,
    PicardAlgebra,
    actual_schedule_picard_certificate,
    iterate_natural_picard,
)
from openai_ns_reconstruction.natural_scale_selection_wide import (
    select_natural_scale_wide,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData


def _schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_actual_schedule_closes_the_two_axis_contraction_scalar_gates() -> None:
    cert = actual_schedule_picard_certificate(_schedule_data(), 0.05)

    # These are exactly the two scalar hypotheses consumed by
    # AxisContraction.exists_fixedPoint_of_controlled after s=1/(2*Lambda).
    assert cert.one_step_radius_upper <= Decimal(1)
    assert cert.contraction_factor_upper <= Decimal("0.5")
    assert cert.Lambda.is_finite()
    assert cert.inverse_two_lambda_upper > 0

    # The wide certificate stays meaningful even though the landed conservative
    # constants are far outside ordinary binary64 range.
    assert cert.remainder_bound_upper.is_finite()
    assert cert.remainder_lipschitz_upper.is_finite()


def test_contraction_certificate_rejects_a_scale_below_its_pinned_threshold() -> None:
    scale = select_natural_scale_wide(
        remainder_bound=Decimal(2),
        remainder_lipschitz=Decimal(3),
        phase_real_part_sup=Decimal("0.1"),
    )
    broken = type(scale)(
        remainder_bound=scale.remainder_bound,
        remainder_lipschitz=scale.remainder_lipschitz,
        phase_real_part_sup=scale.phase_real_part_sup,
        contraction=scale.contraction,
        stability=scale.stability,
        Lambda=scale.contraction - Decimal(1),
        C=scale.C,
    )

    with pytest.raises(ValueError, match="dominate the pinned contraction threshold"):
        NaturalPicardContractionCertificate.from_scale(broken)


def test_tail_bound_matches_an_independent_geometric_sum() -> None:
    scale = select_natural_scale_wide(
        remainder_bound=Decimal(2),
        remainder_lipschitz=Decimal(3),
        phase_real_part_sup=Decimal(0),
    )
    cert = NaturalPicardContractionCertificate.from_scale(scale)

    n = 4
    first = Decimal("0.125")
    bound = cert.tail_bound_from_first_step(first, n)

    with localcontext() as ctx:
        ctx.prec = 90
        q = cert.contraction_factor_upper
        oracle = (q**n) / (Decimal(1) - q) * first
    assert bound >= oracle
    assert bound - oracle < Decimal("1e-80")


def test_execution_adapter_uses_the_pinned_picard_map_shape() -> None:
    # This Decimal scalar fixture checks only the generic map plumbing.  It is
    # deliberately not presented as a coefficient-space or paper-exact field.
    scale = select_natural_scale_wide(
        remainder_bound=Decimal(1),
        remainder_lipschitz=Decimal(1),
        phase_real_part_sup=Decimal(0),
    )
    cert = NaturalPicardContractionCertificate.from_scale(scale)
    reference = Decimal("0.25")

    def remainder(x: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = 90
            return Decimal(1) + x

    def scale_inverse_two_lambda(value: Decimal, Lambda: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = 90
            return value / (Decimal(2) * Lambda)

    def add(x: Decimal, y: Decimal) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = 90
            return x + y

    states = iterate_natural_picard(
        reference,
        certificate=cert,
        remainder=remainder,
        algebra=PicardAlgebra(
            add=add,
            scale_inverse_two_lambda=scale_inverse_two_lambda,
        ),
        iterations=3,
    )

    with localcontext() as ctx:
        ctx.prec = 90
        expected = [reference]
        for _ in range(3):
            expected.append(
                reference
                + (Decimal(1) + expected[-1]) / (Decimal(2) * cert.Lambda)
            )
    assert states == tuple(expected)


def test_iteration_and_tail_interfaces_fail_closed_on_bad_counts() -> None:
    scale = select_natural_scale_wide(
        remainder_bound=Decimal(1),
        remainder_lipschitz=Decimal(1),
        phase_real_part_sup=Decimal(0),
    )
    cert = NaturalPicardContractionCertificate.from_scale(scale)

    with pytest.raises(ValueError, match="nonnegative integer"):
        cert.tail_multiplier_upper(-1)
    with pytest.raises(ValueError, match="nonnegative integer"):
        iterate_natural_picard(
            Decimal(0),
            certificate=cert,
            remainder=lambda x: x,
            algebra=PicardAlgebra(
                add=lambda x, y: x + y,
                scale_inverse_two_lambda=lambda value, Lambda: value / (2 * Lambda),
            ),
            iterations=-1,
        )
