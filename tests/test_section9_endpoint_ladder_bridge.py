import pytest

from openai_ns_reconstruction.section9_endpoint_bridge import (
    Section9UniformEndpointDerivativeWitness,
)
from openai_ns_reconstruction.section9_endpoint_ladder_bridge import (
    admit_section9_uniform_endpoint_ladder,
)
from openai_ns_reconstruction.section9_flat_remainder import (
    Section9FlatRemainderPowerWitness,
    certify_flat_remainder_power_ladder,
)
from openai_ns_reconstruction.section9_stage_certificate import CertifiedBoundDatum


def _endpoint_witness(degree, *, stage=14, spatial_window=3, valid_from=None):
    if valid_from is None:
        valid_from = 0.75 + 0.02 * degree
    return Section9UniformEndpointDerivativeWitness(
        endpoint_derivative_degree=degree,
        section9_derivative_order=degree + 1,
        stage=stage,
        spatial_window=spatial_window,
        coefficient_bound=CertifiedBoundDatum(
            upper_bound=2.0 + degree,
            kind="paper-derived",
            provenance=f"uniform late-time Section 9 derivative theorem for degree {degree}",
        ),
        singularity_exponent=0.20 + 0.05 * degree,
        valid_from=valid_from,
    )


def test_uniform_endpoint_ladder_normalizes_degrees_and_preserves_fail_closed_flags():
    record = admit_section9_uniform_endpoint_ladder(
        [_endpoint_witness(2), _endpoint_witness(0), _endpoint_witness(1)],
        max_endpoint_degree=2,
    )

    assert record.formal_bridge_ready
    assert record.endpoint_degrees == (0, 1, 2)
    assert tuple(row.section9_derivative_order for row in record.bridge_records) == (1, 2, 3)
    assert record.stage == 14
    assert record.spatial_window == 3
    assert record.endpoint_ladder.max_degree == 2
    assert record.endpoint_ladder.spatial_window == 3
    assert record.endpoint_ladder.common_valid_from == pytest.approx(0.79)
    assert record.section9_flat_remainder_ladder_sufficient is False
    assert record.source_theorems_machine_verified is False
    assert record.actual_section9_sequence_verified is False
    assert record.endpoint_limits_constructed is False
    assert record.paper_exact_velocity_available is False


@pytest.mark.parametrize(
    "rows,max_degree,pattern",
    [
        (lambda: [_endpoint_witness(0), _endpoint_witness(2)], 2, "cover exactly"),
        (lambda: [_endpoint_witness(0), _endpoint_witness(1), _endpoint_witness(1)], 1, "duplicate"),
        (lambda: [_endpoint_witness(0), _endpoint_witness(1), _endpoint_witness(2)], 1, "cover exactly"),
    ],
)
def test_uniform_endpoint_ladder_rejects_missing_duplicate_or_extra_degrees(rows, max_degree, pattern):
    with pytest.raises(ValueError, match=pattern):
        admit_section9_uniform_endpoint_ladder(rows(), max_endpoint_degree=max_degree)


def test_uniform_endpoint_ladder_rejects_mixed_stage_or_spatial_window():
    with pytest.raises(ValueError, match="one Section 9 stage"):
        admit_section9_uniform_endpoint_ladder(
            [_endpoint_witness(0, stage=14), _endpoint_witness(1, stage=15)],
            max_endpoint_degree=1,
        )

    with pytest.raises(ValueError, match="one spatial window"):
        admit_section9_uniform_endpoint_ladder(
            [
                _endpoint_witness(0, spatial_window=3),
                _endpoint_witness(1, spatial_window=4),
            ],
            max_endpoint_degree=1,
        )


def test_finite_q_flat_ladder_cannot_bypass_uniform_physical_time_witnesses():
    q_rows = [
        Section9FlatRemainderPowerWitness(
            stage=14,
            derivative_order=1,
            power=power,
            q_upper=0.125,
            constant_bound=CertifiedBoundDatum(
                upper_bound=3.0 + power,
                kind="paper-derived",
                provenance=f"Lemma 9.8 q-flat theorem for N={power}",
            ),
        )
        for power in (0, 1)
    ]
    q_ladder = certify_flat_remainder_power_ladder(q_rows, max_power=1)
    assert q_ladder.formal_ladder_ready

    with pytest.raises(TypeError, match="q-flat"):
        admit_section9_uniform_endpoint_ladder(q_ladder, max_endpoint_degree=1)
