from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import (
    FiniteHierarchyBoundCertificate,
    HierarchyBoundDependency,
    HierarchyJetBound,
    build_slow_borel_schedule_from_hierarchy_certificate,
)
from openai_ns_reconstruction.background_prefix_extension import (
    certify_finite_coherent_slow_borel_residual_chain,
    certify_successive_slow_borel_residual_extension,
)
from openai_ns_reconstruction.background_recurrence_cancellation import (
    ExactRetainedRecurrenceIdentity,
    FiniteRetainedCancellationCertificate,
    FormalAtomTerm,
    RetainedRecurrenceDependency,
    certify_full_residual_tail_majorant,
)
from openai_ns_reconstruction.background_target_decay import (
    certify_finite_physical_power_ladder,
)
from openai_ns_reconstruction.background_truncation_residual import slow_order_exact


H = 0.01
HIERARCHY = "section5-target-fixture-hierarchy"
REVISION = "target-fixture-revision"
STATE = "target-fixture-coefficient-state"


def _bound_certificate(max_order: int):
    entries = []
    for j in range(1, max_order + 1):
        for m in range(j + 3):
            entries.append(
                HierarchyJetBound(
                    order=j,
                    derivative_order=m,
                    bound=float(2 + j + m),
                    hierarchy_id=HIERARCHY,
                    source_revision=REVISION,
                    dependencies=(
                        HierarchyBoundDependency(
                            coefficient_order=j,
                            artifact=f"target-fixture-coefficient-{j}",
                            provider="target-fixture-bound-provider",
                        ),
                    ),
                )
            )
    return FiniteHierarchyBoundCertificate(
        hierarchy_id=HIERARCHY,
        source_revision=REVISION,
        max_order=max_order,
        entries=tuple(entries),
    )


def _cutoff(max_order: int):
    return build_slow_borel_schedule_from_hierarchy_certificate(
        H,
        _bound_certificate(max_order),
        initial_lower_bound=2,
    )


def _identity(order: int):
    atom = f"target-row-{order}"
    return ExactRetainedRecurrenceIdentity(
        order=order,
        hierarchy_id=HIERARCHY,
        source_revision=REVISION,
        coefficient_state_id=STATE,
        identity_source="target fixture exact algebra",
        theorem_name=f"target_fixture_recurrence_{order}",
        dependencies=(
            RetainedRecurrenceDependency(
                coefficient_order=order,
                artifact=f"target-fixture-coefficient-{order}",
                provider="target-fixture-recurrence-provider",
            ),
        ),
        linear_terms=(FormalAtomTerm(atom, Fraction(1, 1)),),
        pair_terms=(FormalAtomTerm(atom, Fraction(-1, 1)),),
        previous_shifted_terms=(),
    )


def _residual(max_order: int, *, q: float = 0.1):
    cancellations = FiniteRetainedCancellationCertificate(
        hierarchy_id=HIERARCHY,
        source_revision=REVISION,
        coefficient_state_id=STATE,
        max_order=max_order,
        identities=tuple(_identity(n) for n in range(max_order + 1)),
    )
    size = max_order + 1
    pair = [[0.0 for _ in range(size)] for _ in range(size)]
    shifted = [0.0 for _ in range(size)]
    pair[max_order][max_order] = 1.0
    shifted[max_order] = 1.0
    return certify_full_residual_tail_majorant(
        cancellations,
        q=q,
        h=H,
        base_power=1.0,
        pair=pair,
        shifted=shifted,
    )


def _extension(previous_order: int, *, q: float = 0.1):
    return certify_successive_slow_borel_residual_extension(
        _cutoff(previous_order),
        _cutoff(previous_order + 1),
        _residual(previous_order, q=q),
        _residual(previous_order + 1, q=q),
    )


def _chain(*, q: float = 0.1):
    return certify_finite_coherent_slow_borel_residual_chain(
        (_extension(1, q=q), _extension(2, q=q))
    )


def test_finite_power_ladder_binds_exact_physical_exponents_and_prefactors():
    chain = _chain()
    targets = (Fraction(1, 1), Fraction(51, 50), Fraction(26, 25))

    certificate = certify_finite_physical_power_ladder(chain, targets)

    endpoints = certificate.endpoints
    expected_exponents = tuple(
        Fraction.from_float(endpoint.tail.base_power)
        + slow_order_exact(endpoint.tail.h, endpoint.first_omitted_order)
        for endpoint in endpoints
    )
    assert certificate.physical_exponents_exact == expected_exponents
    assert certificate.exponent_margins_exact == tuple(
        exponent - target for exponent, target in zip(expected_exponents, targets)
    )
    assert all(value == value for value in certificate.normalized_log_majorants)
    assert certificate.finite_prefix_only is True
    assert certificate.infinite_coherent_family is False
    assert certificate.all_jets_flat is False
    assert certificate.super_algebraic is False
    assert certificate.uniform_in_q is False
    assert certificate.paper_exact is False


def test_finite_power_ladder_rejects_target_beyond_first_omitted_exponent():
    chain = _chain()
    with pytest.raises(ValueError, match="does not attain requested physical q-power"):
        certify_finite_physical_power_ladder(
            chain,
            (Fraction(1, 1), Fraction(2, 1), Fraction(3, 1)),
        )


def test_finite_power_ladder_rejects_float_or_nonincreasing_targets():
    chain = _chain()
    with pytest.raises(TypeError, match="floats are forbidden"):
        certify_finite_physical_power_ladder(
            chain,
            (1.0, Fraction(51, 50), Fraction(26, 25)),
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        certify_finite_physical_power_ladder(
            chain,
            (Fraction(1, 1), Fraction(1, 1), Fraction(26, 25)),
        )


def test_finite_power_ladder_rejects_q_equal_one_as_nondecaying_regime():
    chain = _chain(q=1.0)
    with pytest.raises(ValueError, match="0 < q < 1"):
        certify_finite_physical_power_ladder(
            chain,
            (Fraction(1, 1), Fraction(51, 50), Fraction(26, 25)),
        )
