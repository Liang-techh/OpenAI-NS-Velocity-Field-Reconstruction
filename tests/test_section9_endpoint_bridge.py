from fractions import Fraction

import pytest

from openai_ns_reconstruction.section9_endpoint_bridge import (
    Section9UniformEndpointDerivativeWitness,
    admit_section9_uniform_endpoint_witness,
)
from openai_ns_reconstruction.section9_stage_certificate import (
    CertifiedBoundDatum,
    check_eq_9_18_pointwise,
)


def _uniform_coefficient(value: float = 2.0) -> CertifiedBoundDatum:
    return CertifiedBoundDatum(
        value,
        "rigorous-external",
        "uniform interval theorem on K_3 x [3/4,1)",
    )


def test_uniform_order_n_plus_one_witness_enters_official_endpoint_plateau():
    witness = Section9UniformEndpointDerivativeWitness(
        endpoint_derivative_degree=2,
        section9_derivative_order=3,
        stage=7100,
        spatial_window=3,
        coefficient_bound=_uniform_coefficient(2.0),
        singularity_exponent=0.5,
    )
    record = admit_section9_uniform_endpoint_witness(witness)

    assert record.formal_bridge_ready is True
    assert record.majorant.derivative_degree == 2
    assert record.majorant.spatial_window == 3
    assert record.localization_transfer.certified is True
    assert record.localization_transfer.correction_coefficients(0.8) == (1.0, 0.0, 0.0)

    # Independent antiderivative oracle:
    # C/(1-alpha) * (1-t)^(1-alpha)
    # = 2/(1/2) * sqrt(1/16) = 1 at t=15/16.
    assert record.majorant.endpoint_tail_bound(15.0 / 16.0) == pytest.approx(1.0)
    assert record.source_theorem_machine_verified is False
    assert record.actual_section9_sequence_verified is False
    assert record.endpoint_limit_constructed is False
    assert record.paper_exact_velocity_available is False


def test_bridge_rejects_derivative_order_off_by_one_and_transition_window():
    with pytest.raises(ValueError, match="degree \\+ 1"):
        Section9UniformEndpointDerivativeWitness(
            endpoint_derivative_degree=2,
            section9_derivative_order=2,
            stage=10,
            spatial_window=0,
            coefficient_bound=_uniform_coefficient(),
            singularity_exponent=0.25,
        )

    with pytest.raises(ValueError, match="3/4"):
        Section9UniformEndpointDerivativeWitness(
            endpoint_derivative_degree=0,
            section9_derivative_order=1,
            stage=10,
            spatial_window=0,
            coefficient_bound=_uniform_coefficient(),
            singularity_exponent=0.25,
            valid_from=0.70,
        )


def test_bridge_rejects_nonintegrable_or_non_section10_endpoint_witnesses():
    for exponent in (-0.1, 1.0, float("nan")):
        with pytest.raises(ValueError):
            Section9UniformEndpointDerivativeWitness(
                endpoint_derivative_degree=0,
                section9_derivative_order=1,
                stage=10,
                spatial_window=0,
                coefficient_bound=_uniform_coefficient(),
                singularity_exponent=exponent,
            )

    with pytest.raises(ValueError, match="endpoint=1"):
        Section9UniformEndpointDerivativeWitness(
            endpoint_derivative_degree=0,
            section9_derivative_order=1,
            stage=10,
            spatial_window=0,
            coefficient_bound=_uniform_coefficient(),
            singularity_exponent=0.25,
            endpoint=2.0,
        )


def test_pointwise_eq_918_certificate_cannot_be_promoted_to_uniform_endpoint_majorant():
    pointwise = check_eq_9_18_pointwise(
        stage=6998,
        derivative_order=1,
        h=Fraction(1, 200),
        q=0.25,
        derivative_loss=Fraction(3, 2),
        constant=2.0,
        log_power=1,
        residual_bound=CertifiedBoundDatum(
            1.0e-4, "certified-numerical", "one-q interval enclosure"
        ),
        flat_remainder_bound=CertifiedBoundDatum(
            1.0e-6, "certified-numerical", "one-q flat-remainder enclosure"
        ),
    )
    assert pointwise.arithmetic_check_passed is True

    with pytest.raises(TypeError, match="pointwise"):
        admit_section9_uniform_endpoint_witness(pointwise)  # type: ignore[arg-type]


def test_uniform_witness_requires_auditable_bound_provenance():
    with pytest.raises(ValueError, match="kind"):
        CertifiedBoundDatum(2.0, "sampled", "grid maximum")
    with pytest.raises(ValueError, match="provenance"):
        CertifiedBoundDatum(2.0, "paper-derived", "   ")
