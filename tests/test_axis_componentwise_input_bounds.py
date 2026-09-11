import pytest

from openai_ns_reconstruction.axis_componentwise_input_bounds import (
    componentwise_analytic_input_norm_certificate,
)


def _field_bounds(*, chi: float = 2.0, gradient: float = 1.0e12):
    return {
        "one": 1.0,
        "eta": 1.2,
        "d": 2.5,
        "inverseL": 1.1,
        "uStar": 4.5,
        "uStarEta": 4.0,
        "wStar": 8.0,
        "hStar": 3.0,
        "zStar": 20.0,
        "chi": chi,
        "gradient": gradient,
    }


def test_componentwise_ledger_keeps_gradient_out_of_chi_resolvent_bound() -> None:
    cert = componentwise_analytic_input_norm_certificate(
        neighborhood_radius=0.2,
        field_value_sup_upper=_field_bounds(),
    )

    assert cert.epsilon == 0.1
    assert cert.radius_loss_half == 12.0
    assert cert.amplitude_norm_upper == 12.0

    # Positive products are rounded upward, so compare by inequalities rather
    # than assuming the exact binary64 representation of 12*B_k.
    assert cert.field_norm_upper["chi"] >= 24.0
    assert cert.field_norm_upper["gradient"] >= 12.0e12
    assert cert.chi_norm_upper == cert.field_norm_upper["chi"]
    assert cert.coefficient_family_norm_upper == cert.field_norm_upper["gradient"]

    # The AxisResolvent parameter depends on chi specifically, not on the
    # unrelated maximum family bound.  This is the entire point of the refined
    # ledger and is independent of the chosen large gradient regression value.
    K = cert.resolvent_majorant.majorant_parameter_upper
    assert K >= 2560.0 * 24.0
    assert K < 100_000.0
    assert K < 2560.0 * cert.coefficient_family_norm_upper
    assert cert.resolvent_majorant.series_upper > 1.0


def test_componentwise_axis_data_uses_each_matching_fixed_field_bound() -> None:
    cert = componentwise_analytic_input_norm_certificate(
        neighborhood_radius=0.2,
        field_value_sup_upper=_field_bounds(chi=1.5, gradient=7.0),
    )
    data = cert.axis_data_norm_bounds(h=0.01)
    n = cert.field_norm_upper

    assert data.A == pytest.approx(0.51)
    assert data.D == pytest.approx(0.49)
    assert data.one == n["one"]
    assert data.eta == n["eta"]
    assert data.d == n["d"]
    assert data.inverse_l == n["inverseL"]
    assert data.u_star == n["uStar"]
    assert data.u_star_eta == n["uStarEta"]
    assert data.w_star == n["wStar"]
    assert data.h_star == n["hStar"]
    assert data.normalized_gradient == n["gradient"]
    assert data.z_star == n["zStar"]


def test_componentwise_ledger_fails_closed_on_missing_or_invalid_fields() -> None:
    missing = _field_bounds()
    missing.pop("chi")
    with pytest.raises(ValueError, match="pinned fixed fields"):
        componentwise_analytic_input_norm_certificate(
            neighborhood_radius=0.2,
            field_value_sup_upper=missing,
        )

    invalid = _field_bounds()
    invalid["chi"] = 0.0
    with pytest.raises(ValueError, match="field value bound chi"):
        componentwise_analytic_input_norm_certificate(
            neighborhood_radius=0.2,
            field_value_sup_upper=invalid,
        )
