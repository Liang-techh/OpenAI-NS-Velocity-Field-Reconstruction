from fractions import Fraction

import pytest

from openai_ns_reconstruction.background_hierarchy_cutoff_bounds import (
    FiniteHierarchyBoundCertificate,
    HierarchyBoundDependency,
    HierarchyJetBound,
    build_slow_borel_schedule_from_hierarchy_certificate,
)
from openai_ns_reconstruction.background_prefix_extension import (
    certify_successive_slow_borel_residual_extension,
)
from openai_ns_reconstruction.background_recurrence_cancellation import (
    ExactRetainedRecurrenceIdentity,
    FiniteRetainedCancellationCertificate,
    FormalAtomTerm,
    RetainedRecurrenceDependency,
    certify_full_residual_tail_majorant,
)


H = 0.01
HIERARCHY = "section5-fixture-hierarchy"
REVISION = "fixture-revision"
STATE = "fixture-coefficient-state"


def _bound_certificate(max_order: int, *, first_bound_override: float | None = None):
    entries = []
    for j in range(1, max_order + 1):
        for m in range(j + 3):
            bound = float(2 + j + m)
            if j == 1 and m == 0 and first_bound_override is not None:
                bound = float(first_bound_override)
            entries.append(
                HierarchyJetBound(
                    order=j,
                    derivative_order=m,
                    bound=bound,
                    hierarchy_id=HIERARCHY,
                    source_revision=REVISION,
                    dependencies=(
                        HierarchyBoundDependency(
                            coefficient_order=j,
                            artifact=f"fixture-coefficient-{j}",
                            provider="fixture-analytic-bound-provider",
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


def _cutoff(max_order: int, *, first_bound_override: float | None = None):
    return build_slow_borel_schedule_from_hierarchy_certificate(
        H,
        _bound_certificate(max_order, first_bound_override=first_bound_override),
        initial_lower_bound=2,
    )


def _identity(order: int, *, theorem_name: str | None = None):
    atom = f"fixture-row-{order}"
    return ExactRetainedRecurrenceIdentity(
        order=order,
        hierarchy_id=HIERARCHY,
        source_revision=REVISION,
        coefficient_state_id=STATE,
        identity_source="fixture exact algebra",
        theorem_name=theorem_name or f"fixture_recurrence_{order}",
        dependencies=(
            RetainedRecurrenceDependency(
                coefficient_order=order,
                artifact=f"fixture-coefficient-{order}",
                provider="fixture-recurrence-provider",
            ),
        ),
        linear_terms=(FormalAtomTerm(atom, Fraction(1, 1)),),
        pair_terms=(FormalAtomTerm(atom, Fraction(-1, 1)),),
        previous_shifted_terms=(),
    )


def _cancellations(max_order: int, *, replacement_order1=None):
    identities = [_identity(n) for n in range(max_order + 1)]
    if replacement_order1 is not None and max_order >= 1:
        identities[1] = replacement_order1
    return FiniteRetainedCancellationCertificate(
        hierarchy_id=HIERARCHY,
        source_revision=REVISION,
        coefficient_state_id=STATE,
        max_order=max_order,
        identities=tuple(identities),
    )


def _residual(max_order: int, *, q: float = 0.1, large_new_prefactor: bool = False):
    cancellations = _cancellations(max_order)
    size = max_order + 1
    pair = [[0.0 for _ in range(size)] for _ in range(size)]
    shifted = [0.0 for _ in range(size)]
    if max_order == 1:
        pair[1][1] = 3.0
        shifted[1] = 1.0
    elif max_order == 2:
        if large_new_prefactor:
            pair[2][2] = 100.0
            shifted[2] = 100.0
        else:
            pair[2][2] = 1.0
            shifted[2] = 1.0
    else:
        raise ValueError("fixture only defines orders one and two")
    return certify_full_residual_tail_majorant(
        cancellations,
        q=q,
        h=H,
        base_power=1.0,
        pair=pair,
        shifted=shifted,
    )


def test_successive_extension_binds_cutoff_recurrence_and_tail_improvement():
    previous_cutoff = _cutoff(1)
    extended_cutoff = _cutoff(2)
    previous_residual = _residual(1)
    extended_residual = _residual(2)

    witness = certify_successive_slow_borel_residual_extension(
        previous_cutoff,
        extended_cutoff,
        previous_residual,
        extended_residual,
    )

    assert witness.previous_order == 1
    assert witness.extended_order == 2
    assert witness.previous_first_omitted_order == 2
    assert witness.extended_first_omitted_order == 3
    assert witness.exponent_gain_exact == 2 * Fraction.from_float(H)
    assert witness.new_recursive_scale == extended_cutoff.schedule.scales[-1]
    assert witness.majorant_nonincreasing
    assert witness.majorant_strictly_improves
    assert witness.paper_exact is False


def test_extension_rejects_changed_old_hierarchy_bound():
    with pytest.raises(ValueError, match="changes a previously certified bound"):
        certify_successive_slow_borel_residual_extension(
            _cutoff(1),
            _cutoff(2, first_bound_override=9.0),
            _residual(1),
            _residual(2),
        )


def test_extension_rejects_changed_old_exact_recurrence_identity():
    previous_residual = _residual(1)
    changed = _identity(1, theorem_name="different_exact_theorem_name")
    new_cancellations = _cancellations(2, replacement_order1=changed)
    size = 3
    pair = [[0.0 for _ in range(size)] for _ in range(size)]
    pair[2][2] = 1.0
    shifted = [0.0, 0.0, 1.0]
    extended_residual = certify_full_residual_tail_majorant(
        new_cancellations,
        q=0.1,
        h=H,
        base_power=1.0,
        pair=pair,
        shifted=shifted,
    )

    with pytest.raises(ValueError, match="changes a retained identity"):
        certify_successive_slow_borel_residual_extension(
            _cutoff(1),
            _cutoff(2),
            previous_residual,
            extended_residual,
        )


def test_extension_rejects_different_q():
    with pytest.raises(ValueError, match="same q"):
        certify_successive_slow_borel_residual_extension(
            _cutoff(1),
            _cutoff(2),
            _residual(1, q=0.1),
            _residual(2, q=0.2),
        )


def test_extension_rejects_exponent_gain_without_majorant_improvement():
    with pytest.raises(ValueError, match="majorant is larger"):
        certify_successive_slow_borel_residual_extension(
            _cutoff(1),
            _cutoff(2),
            _residual(1),
            _residual(2, large_new_prefactor=True),
        )
