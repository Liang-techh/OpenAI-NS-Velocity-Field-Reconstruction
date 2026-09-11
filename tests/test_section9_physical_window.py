from fractions import Fraction

import numpy as np
import pytest

from openai_ns_reconstruction.coordinates import solve_q
from openai_ns_reconstruction.section9_local_sum_admission import Section9LocalSumHypothesis
from openai_ns_reconstruction.section9_physical_window import (
    certify_section10_physical_slab_local_finiteness,
)


def _hypothesis(q_big=Fraction(1, 2)):
    return Section9LocalSumHypothesis(
        q_big=q_big,
        first_scale=8,
        kind="paper-derived",
        provenance="Lemma 9.7 common-domain + Lemma 5.4 infinite doubling hypotheses",
    )


def test_physical_slab_maps_to_exact_q_strip_and_local_sum_bound():
    cert = certify_section10_physical_slab_local_finiteness(
        _hypothesis(),
        t_lower=Fraction(3, 4),
        t_upper=Fraction(15, 16),
    )

    assert cert.q_lower_bound == Fraction(1, 16)
    assert cert.q_upper_bound == Fraction(5, 16)
    assert cert.local_sum.q_lower == cert.q_lower_bound
    assert cert.local_sum.q_upper == cert.q_upper_bound
    assert cert.max_potentially_active_stage == 1
    assert cert.first_guaranteed_inactive_stage == 2
    assert cert.support_geometry_verified
    assert cert.similarity_bound_verified
    assert not cert.endpoint_covered
    assert not cert.actual_correction_fields_verified
    assert not cert.eq_9_21_sum_constructed
    assert not cert.endpoint_uniform_residual_majorants_verified
    assert not cert.paper_exact_velocity_available


def test_independent_similarity_solver_stays_inside_analytic_bounds():
    cert = certify_section10_physical_slab_local_finiteness(
        _hypothesis(),
        t_lower=Fraction(3, 4),
        t_upper=Fraction(15, 16),
    )

    q_lo = float(cert.q_lower_bound)
    q_hi = float(cert.q_upper_bound)
    for h in (0.005, 0.1, 0.2):
        for t in np.linspace(0.75, 15.0 / 16.0, 7):
            for z in np.linspace(-0.25, 0.25, 9):
                q = solve_q(float(z), float(t), h)
                assert q >= q_lo * (1.0 - 2e-12)
                assert q <= q_hi * (1.0 + 2e-12)


def test_later_slab_tightens_q_upper_bound_without_sampling_q_big():
    cert = certify_section10_physical_slab_local_finiteness(
        _hypothesis(Fraction(1, 4)),
        t_lower=Fraction(7, 8),
        t_upper=Fraction(15, 16),
    )
    assert cert.q_lower_bound == Fraction(1, 16)
    assert cert.q_upper_bound == Fraction(3, 16)
    assert cert.q_upper_bound < cert.local_sum.q_big


def test_bridge_fails_closed_when_common_q_domain_is_too_small():
    with pytest.raises(ValueError, match="not inside q_big"):
        certify_section10_physical_slab_local_finiteness(
            _hypothesis(Fraction(1, 4)),
            t_lower=Fraction(3, 4),
            t_upper=Fraction(15, 16),
        )


def test_bridge_does_not_promote_a_finite_slab_to_the_endpoint():
    with pytest.raises(ValueError, match="strictly before t=1"):
        certify_section10_physical_slab_local_finiteness(
            _hypothesis(),
            t_lower=Fraction(3, 4),
            t_upper=1,
        )


def test_bridge_requires_official_late_plateau_and_ordered_slab():
    with pytest.raises(ValueError, match="official t>=3/4 plateau"):
        certify_section10_physical_slab_local_finiteness(
            _hypothesis(),
            t_lower=Fraction(2, 3),
            t_upper=Fraction(7, 8),
        )
    with pytest.raises(ValueError, match="t_lower < t_upper"):
        certify_section10_physical_slab_local_finiteness(
            _hypothesis(),
            t_lower=Fraction(7, 8),
            t_upper=Fraction(7, 8),
        )
