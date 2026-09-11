import numpy as np
import pytest
from scipy.integrate import quad

from openai_ns_reconstruction.mean_rank_update import (
    PAPER_EXACT,
    STATUS,
    MeanRankUpdatePlan,
    cell_intervals,
    inner_support,
)


def _integral_over_supports(function, supports):
    return sum(
        quad(function, left, right, epsabs=1e-10, epsrel=1e-10, limit=200)[0]
        for left, right in supports
    )


def test_pinned_cell_geometry_and_debt_scaling():
    # For a=1,b=8 the 3-cell step is exactly 1.  This hand-computed geometry
    # is independent of the production moment matrix and bump evaluation.
    assert cell_intervals(1.0, 8.0, 3) == (
        (2.0, 3.0),
        (4.0, 5.0),
        (6.0, 7.0),
    )
    assert inner_support((2.0, 3.0)) == (2.25, 2.75)

    plan = MeanRankUpdatePlan(
        lam=0.5,
        C=2.0,
        a=1.0,
        b=8.0,
        ell=3.0,
        U=2.0,
        debt=(8.0, 108.0, 72.0),
    )
    # normalizeDebt = (d0/U^2, d1/(ell^3 U^2), d2/(ell^2 U^2)).
    np.testing.assert_allclose(plan.normalized_debt, [2.0, 1.0, 2.0])
    np.testing.assert_allclose(plan.angular_targets, [0.0, -0.5, 1.0])
    np.testing.assert_allclose(plan.axial_targets, [0.0, -0.5])
    np.testing.assert_allclose(
        plan.rows_from_normalized_moments(plan.angular_targets, plan.axial_targets),
        [0.0, 0.0, -8.0, -108.0, -72.0],
        rtol=0.0,
        atol=1e-13,
    )


def test_compact_mean_correction_matches_independent_physical_five_rows():
    """Verify the solved correction with SciPy quadrature, not its solve matrix."""
    plan = MeanRankUpdatePlan(
        lam=0.5,
        C=2.0,
        a=1.0,
        b=8.0,
        ell=3.0,
        U=2.0,
        debt=(8.0, 108.0, 72.0),
        quadrature_points=256,
    )
    correction = plan.solve()

    dv = correction.angular_increment
    ga = correction.desired_axial_increment
    V = correction.background
    angular_supports = correction.physical_angular_supports
    axial_supports = correction.physical_axial_supports

    rows = np.array(
        [
            _integral_over_supports(lambda r: r**2 * dv(r), angular_supports),
            _integral_over_supports(lambda r: r * ga(r), axial_supports),
            _integral_over_supports(lambda r: (2.0 * V(r) / r) * dv(r), angular_supports),
            _integral_over_supports(lambda r: r**2 * V(r) * ga(r), axial_supports),
            _integral_over_supports(lambda r: -r * V(r) * dv(r), angular_supports),
        ]
    )
    np.testing.assert_allclose(rows, plan.physical_row_targets, rtol=2e-10, atol=2e-9)

    # Exact compact-support behavior of the executable representative.
    for left, right in angular_supports:
        assert correction.angular_increment(left) == 0.0
        assert correction.angular_increment(right) == 0.0
    for left, right in axial_supports:
        assert correction.desired_axial_increment(left) == 0.0
        assert correction.desired_axial_increment(right) == 0.0


def test_mean_rank_update_stays_formal_structure_and_fails_closed():
    assert STATUS == "formal-structure"
    assert PAPER_EXACT is False
    plan = MeanRankUpdatePlan(0.5, 2.0, 1.0, 8.0, 3.0, 2.0, (1.0, 2.0, 3.0))
    assert plan.paper_exact is False
    assert plan.solve().paper_exact is False

    with pytest.raises(ValueError, match="lam must be positive"):
        MeanRankUpdatePlan(0.0, 2.0, 1.0, 8.0, 3.0, 2.0, (1.0, 2.0, 3.0))
    with pytest.raises(ValueError, match="C must be nonzero"):
        MeanRankUpdatePlan(0.5, 0.0, 1.0, 8.0, 3.0, 2.0, (1.0, 2.0, 3.0))
    with pytest.raises(ValueError, match="ell must be positive"):
        MeanRankUpdatePlan(0.5, 2.0, 1.0, 8.0, 0.0, 2.0, (1.0, 2.0, 3.0))
    with pytest.raises(ValueError, match="U must be nonzero"):
        MeanRankUpdatePlan(0.5, 2.0, 1.0, 8.0, 3.0, 0.0, (1.0, 2.0, 3.0))
