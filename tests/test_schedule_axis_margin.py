import numpy as np
import pytest

from openai_ns_reconstruction.natural_axis import H
from openai_ns_reconstruction.natural_axis_range import (
    NaturalAxisRangeParameters,
    locate_unique_H_root,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.schedule_axis_margin import (
    certify_schedule_low_Z_margin,
    clock_mass_upper,
    release_start_upper,
    schedule_natural_axis_Z,
)
from openai_ns_reconstruction.schedule_axis_pressure import schedule_pressure_mass


def _data(*, P: float = 2.0, h: float = 0.01) -> TailData:
    return TailData(
        OutgoingCoreParameters(P=P, m=1.0, lam=0.05, wait=30.0),
        h=h,
    )


def test_schedule_envelope_bounds_actual_clock_mass_without_tail_quadrature_assumption() -> None:
    data = _data()
    mass_upper, release_upper, rho_upper = clock_mass_upper(data)

    assert release_upper > data.release_start
    assert 0.0 < rho_upper < 1.0
    # At eta=0 the pressure mass is exactly the total clockWeight mass.  The
    # production evaluator is not used to *construct* the analytic upper bound;
    # this is an independent numerical cross-check of the resulting inequality.
    assert schedule_pressure_mass(data, 0.0) < mass_upper


def test_actual_schedule_Z_has_the_pinned_root_separation() -> None:
    data = _data()
    parameters = NaturalAxisRangeParameters(h=data.h, j=0.05)
    root = locate_unique_H_root(parameters)
    z_root = schedule_natural_axis_Z(data, parameters.j, root.eta0)

    assert abs(root.residual) < 1e-12
    # NaturalAxisRange.Z_at_root_lower gives the strict theorem-side j/5 gap
    # once SchedulePressure.pressureData is instantiated with P>=2.
    assert z_root > parameters.j / 5.0


def test_constructed_margin_instantiates_sigma_and_survives_pointwise_crosscheck() -> None:
    data = _data()
    witness = certify_schedule_low_Z_margin(data, 0.05)

    assert witness.margin > 0.0
    assert witness.cutoff.delta == pytest.approx(0.005)
    assert witness.cutoff.sigma > 0.0
    assert witness.cutoff.margin == witness.margin
    assert witness.z_lipschitz_upper > 0.0
    assert witness.root_distance_lower > 0.0

    # Sampling is deliberately only a regression cross-check; the witness above
    # comes from the analytic root-separation/Lipschitz chain, not this grid.
    threshold = witness.parameters.j / 10.0
    for eta in np.linspace(-1.0, 1.0, 33):
        eta = float(eta)
        z = schedule_natural_axis_Z(data, witness.parameters.j, eta)
        if abs(z) <= threshold:
            assert H(witness.parameters.h, witness.parameters.j, eta) ** 2 >= witness.margin


def test_margin_constructor_fails_closed_outside_the_actual_theorem_hypotheses() -> None:
    with pytest.raises(ValueError, match="P >= 2"):
        certify_schedule_low_Z_margin(_data(P=1.5), 0.05)

    with pytest.raises(ValueError, match="1/100"):
        certify_schedule_low_Z_margin(_data(h=0.02), 0.05)

    with pytest.raises(ValueError, match="1/20"):
        certify_schedule_low_Z_margin(_data(), 0.051)

    with pytest.raises(TypeError):
        release_start_upper(object())
