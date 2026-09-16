from fractions import Fraction as F

import pytest

from openai_ns_reconstruction.kokuno_profile_material_transport import (
    FULL_RECONSTRUCTION,
    IMPORTED_PROFILE_EXISTENCE_PROVED_HERE,
    PAPER_EXACT_VELOCITY_AVAILABLE,
    ScalarJet,
    TransportPoint,
    coefficient_residuals,
    composed_material_transport_core,
    material_transport_reduction_residual,
    radial_primitive_value,
    reduced_material_transport_core,
    transport_Hc,
    transport_W,
)


def fixture():
    point = TransportPoint(h=F(1, 8), eta=F(1, 3), X=F(5, 7))
    f = ScalarJet(value=F(7, 5), d_eta=F(-4, 11), d_X=F(13, 17))
    return point, f


def test_exact_public_transport_reduction():
    point, f = fixture()
    U, M, M_eta, b = F(11, 9), F(2, 5), F(-3, 10), F(-2, 3)
    V0 = radial_primitive_value(point, U=U, M=M, M_eta=M_eta)

    assert V0 == F(566, 735)
    assert transport_W(point, M=M, M_eta=M_eta) == F(37, 30)
    assert transport_Hc(point, U=U) == F(785, 648)
    assert composed_material_transport_core(point, b=b, f=f, U=U, V0=V0) == F(19132, 45815)
    assert reduced_material_transport_core(point, b=b, f=f, U=U, M=M, M_eta=M_eta) == F(19132, 45815)
    assert material_transport_reduction_residual(point, b=b, f=f, U=U, M=M, M_eta=M_eta) == 0


def test_each_reduced_coefficient_closes_exactly():
    point, _ = fixture()
    assert coefficient_residuals(
        point,
        b=F(-2, 3),
        U=F(11, 9),
        M=F(2, 5),
        M_eta=F(-3, 10),
    ) == (0, 0, 0)


def test_same_witness_V0_binding_fails_closed_under_tiny_perturbation():
    point, f = fixture()
    U, M, M_eta, b = F(11, 9), F(2, 5), F(-3, 10), F(-2, 3)
    V0 = radial_primitive_value(point, U=U, M=M, M_eta=M_eta)
    eps = F(1, 2**40)
    residual = material_transport_reduction_residual(
        point,
        b=b,
        f=f,
        U=U,
        M=M,
        M_eta=M_eta,
        V0=V0 + eps,
    )
    assert residual == eps * f.d_X
    assert residual == F(13, 18691697672192)
    assert residual != 0


def test_axial_X_chain_term_is_observable():
    point, f = fixture()
    U, M, M_eta, b = F(11, 9), F(2, 5), F(-3, 10), F(-2, 3)
    V0 = radial_primitive_value(point, U=U, M=M, M_eta=M_eta)
    correct = composed_material_transport_core(point, b=b, f=f, U=U, V0=V0)

    time_n = -b * f.value + point.D * point.eta * f.d_eta + point.X * f.d_X
    radial_n = point.L * V0 * f.d_X
    bad_axial_n = U * (2 * b * point.eta * f.value + point.d * f.d_eta)
    bad = (time_n + radial_n + bad_axial_n) / point.L
    defect = bad - correct

    assert defect == (2 * point.eta * point.X * U * f.d_X) / point.L
    assert defect != 0


def test_domain_and_exact_input_guards():
    with pytest.raises(ValueError):
        TransportPoint(h=F(1, 8), eta=F(1, 3), X=0)
    with pytest.raises(ValueError):
        TransportPoint(h=F(1, 2), eta=0, X=1)
    with pytest.raises(TypeError):
        TransportPoint(h=0.125, eta=F(1, 3), X=1)
    with pytest.raises(TypeError):
        ScalarJet(value=1.0, d_eta=0, d_X=0)


def test_truth_boundary_stays_fail_closed():
    assert IMPORTED_PROFILE_EXISTENCE_PROVED_HERE is False
    assert PAPER_EXACT_VELOCITY_AVAILABLE is False
    assert FULL_RECONSTRUCTION is False
