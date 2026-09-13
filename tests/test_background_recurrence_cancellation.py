from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_recurrence_cancellation import (
    ExactRetainedRecurrenceIdentity,
    FiniteRetainedCancellationCertificate,
    FormalAtomTerm,
    RetainedRecurrenceDependency,
    certify_full_residual_tail_majorant,
)


def _dep(order: int) -> RetainedRecurrenceDependency:
    return RetainedRecurrenceDependency(
        coefficient_order=order,
        artifact=f"repaired-coefficient-{order}",
        provider=f"section5.hierarchy[{order}]",
    )


def _identity(order: int, *, state: str = "state-v1") -> ExactRetainedRecurrenceIdentity:
    own = f"row-{order}"
    carry = f"carry-{order}"
    dependencies = (_dep(order),)
    if order:
        dependencies += (_dep(order - 1),)
    return ExactRetainedRecurrenceIdentity(
        order=order,
        hierarchy_id="section5-hierarchy",
        source_revision="source-revision-1",
        coefficient_state_id=state,
        identity_source=f"hierarchy.recurrence[{order}]",
        theorem_name="SlowExpansionResidual.recurrence_truncation_of_zero",
        dependencies=dependencies,
        linear_terms=(
            FormalAtomTerm(own, Fraction(3, 5)),
            FormalAtomTerm(carry, -2),
        ),
        pair_terms=(
            FormalAtomTerm(own, Fraction(2, 5)),
            FormalAtomTerm(carry, 3),
        ),
        previous_shifted_terms=(
            FormalAtomTerm(own, 1),
            FormalAtomTerm(carry, 1),
        ),
    )


def _certificate(order: int = 2) -> FiniteRetainedCancellationCertificate:
    return FiniteRetainedCancellationCertificate(
        hierarchy_id="section5-hierarchy",
        source_revision="source-revision-1",
        coefficient_state_id="state-v1",
        max_order=order,
        identities=tuple(_identity(n) for n in range(order + 1)),
    )


def test_exact_formal_identity_cancels_without_tolerance() -> None:
    identity = _identity(2)
    assert identity.exact_zero is True
    assert identity.residual_terms == ()
    assert identity.order == 2
    assert any(dep.coefficient_order == 2 for dep in identity.dependencies)


def test_tiny_but_nonzero_exact_residual_is_rejected() -> None:
    with pytest.raises(ValueError, match="not an exact zero identity"):
        ExactRetainedRecurrenceIdentity(
            order=1,
            hierarchy_id="section5-hierarchy",
            source_revision="source-revision-1",
            coefficient_state_id="state-v1",
            identity_source="hierarchy.recurrence[1]",
            theorem_name="SlowExpansionResidual.recurrence_truncation_of_zero",
            dependencies=(_dep(1),),
            linear_terms=(FormalAtomTerm("epsilon", Fraction(1, 10**30)),),
            pair_terms=(),
            previous_shifted_terms=(),
        )


def test_float_coefficients_are_forbidden_from_exact_identity() -> None:
    with pytest.raises(TypeError, match="floats are forbidden"):
        FormalAtomTerm("sampled-residual", 1e-30)


def test_dependency_must_bind_own_order_and_never_future_order() -> None:
    with pytest.raises(ValueError, match="own coefficient order"):
        ExactRetainedRecurrenceIdentity(
            order=2,
            hierarchy_id="h",
            source_revision="r",
            coefficient_state_id="s",
            identity_source="i",
            theorem_name="t",
            dependencies=(_dep(1),),
            linear_terms=(),
            pair_terms=(),
            previous_shifted_terms=(),
        )

    with pytest.raises(ValueError, match="future coefficient"):
        ExactRetainedRecurrenceIdentity(
            order=2,
            hierarchy_id="h",
            source_revision="r",
            coefficient_state_id="s",
            identity_source="i",
            theorem_name="t",
            dependencies=(_dep(2), _dep(3)),
            linear_terms=(),
            pair_terms=(),
            previous_shifted_terms=(),
        )


def test_certificate_requires_contiguous_same_state_prefix() -> None:
    with pytest.raises(ValueError, match="contiguous prefix"):
        FiniteRetainedCancellationCertificate(
            hierarchy_id="section5-hierarchy",
            source_revision="source-revision-1",
            coefficient_state_id="state-v1",
            max_order=2,
            identities=(_identity(0), _identity(2)),
        )

    with pytest.raises(ValueError, match="coefficient_state_id"):
        FiniteRetainedCancellationCertificate(
            hierarchy_id="section5-hierarchy",
            source_revision="source-revision-1",
            coefficient_state_id="state-v1",
            max_order=1,
            identities=(_identity(0), _identity(1, state="different-state")),
        )


def test_exact_cancellation_unlocks_first_omitted_order_tail_majorant() -> None:
    certificate = _certificate(2)
    pair = (
        (0.2, -0.1, 0.05),
        (0.07, 0.04, -0.02),
        (-0.03, 0.06, 0.02),
    )
    shifted = (0.11, -0.09, 0.07)
    result = certify_full_residual_tail_majorant(
        certificate,
        q=0.5,
        h=0.125,
        base_power=0.0,
        pair=pair,
        shifted=shifted,
    )

    assert result.first_omitted_order == 3
    assert result.tail.order == certificate.max_order
    assert result.tail.exponent_exact == Fraction(3, 4)
    assert result.tail.coefficient_l1 > 0.0
    assert result.paper_exact is False
    assert result.cancellations.paper_exact is False


def test_tail_gate_rejects_uncertified_objects() -> None:
    pair = ((0.0,),)
    shifted = (0.0,)
    with pytest.raises(TypeError, match="FiniteRetainedCancellationCertificate"):
        certify_full_residual_tail_majorant(
            object(),
            q=0.5,
            h=0.125,
            base_power=0.0,
            pair=pair,
            shifted=shifted,
        )
