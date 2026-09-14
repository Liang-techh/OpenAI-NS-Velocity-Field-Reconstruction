from fractions import Fraction

import pytest

from openai_ns_reconstruction.section8_compact_support_geometry import (
    AffineSupportFamily,
    PINNED_FORMAL_COMMIT,
    Section8CompactSupportGeometry,
    pinned_cell_coordinates,
    pinned_inner_support,
)


def test_section8_pinned_support_geometry_is_definition_exact():
    geometry = Section8CompactSupportGeometry.pinned()

    assert geometry.angular.cells == (
        (Fraction(1, 7), Fraction(2, 7)),
        (Fraction(3, 7), Fraction(4, 7)),
        (Fraction(5, 7), Fraction(6, 7)),
    )
    assert geometry.angular.supports == (
        (Fraction(5, 28), Fraction(1, 4)),
        (Fraction(13, 28), Fraction(15, 28)),
        (Fraction(21, 28), Fraction(23, 28)),
    )
    assert geometry.axial.cells == (
        (Fraction(1, 5), Fraction(2, 5)),
        (Fraction(3, 5), Fraction(4, 5)),
    )
    assert geometry.axial.supports == (
        (Fraction(1, 4), Fraction(7, 20)),
        (Fraction(13, 20), Fraction(3, 4)),
    )

    assert geometry.angular.exact_affine_geometry_verified is True
    assert geometry.angular.pairwise_separated_verified is True
    assert geometry.angular.support_strictly_inside_cells_verified is True
    assert geometry.angular.cells_strictly_inside_patch_verified is True
    assert geometry.axial.exact_affine_geometry_verified is True
    assert geometry.valid_for_any_strict_positive_patch is True
    assert geometry.status == "definition-exact-support-geometry"
    assert geometry.formal_commit == PINNED_FORMAL_COMMIT


def test_support_geometry_matches_independent_affine_formula_replay():
    # Independent direct arithmetic for the first angular cell:
    # cell = (a + 1/7 (b-a), a + 2/7 (b-a)); the repair support is its
    # middle half, hence affine coordinates 5/28 and 7/28.
    cells = pinned_cell_coordinates(3)
    assert cells[0] == (Fraction(1, 7), Fraction(2, 7))
    assert pinned_inner_support(cells[0]) == (Fraction(5, 28), Fraction(7, 28))

    # Exact positive unused gaps are part of the paper construction, not a
    # numerical support-overlap heuristic.
    assert cells[1][0] - cells[0][1] == Fraction(1, 7)
    assert cells[2][0] - cells[1][1] == Fraction(1, 7)


def test_family_rejects_generic_or_rounded_cutoff_geometry():
    exact_cells = pinned_cell_coordinates(3)
    exact_supports = tuple(pinned_inner_support(cell) for cell in exact_cells)

    with pytest.raises(TypeError, match="floats are rejected"):
        AffineSupportFamily(
            count=3,
            cells=((1.0 / 7.0, Fraction(2, 7)), exact_cells[1], exact_cells[2]),
            supports=exact_supports,
        )

    shifted_cells = (
        exact_cells[0],
        (Fraction(22, 49), exact_cells[1][1]),
        exact_cells[2],
    )
    with pytest.raises(ValueError, match="pinned FiveRowRank"):
        AffineSupportFamily(count=3, cells=shifted_cells, supports=exact_supports)

    widened_supports = (
        (Fraction(1, 7), exact_supports[0][1]),
        exact_supports[1],
        exact_supports[2],
    )
    with pytest.raises(ValueError, match="pinned LocalizedMomentRepair"):
        AffineSupportFamily(count=3, cells=exact_cells, supports=widened_supports)


def test_certificate_fails_closed_on_wrong_family_or_provenance():
    angular = AffineSupportFamily.pinned(3)
    axial = AffineSupportFamily.pinned(2)

    with pytest.raises(ValueError, match="exactly three cells"):
        Section8CompactSupportGeometry(
            angular=AffineSupportFamily.pinned(2),
            axial=axial,
        )

    with pytest.raises(ValueError, match="formal_commit"):
        Section8CompactSupportGeometry(
            angular=angular,
            axial=axial,
            formal_commit="caller-supplied-revision",
        )

    geometry = Section8CompactSupportGeometry.pinned()
    assert geometry.actual_patch_endpoints_materialized is False
    assert geometry.actual_bump_values_materialized is False
    assert geometry.exact_moment_identities_machine_replayed is False
    assert geometry.five_row_solve_materialized is False
    assert geometry.actual_wave_defect_consumed is False
    assert geometry.paper_exact_velocity_available is False
