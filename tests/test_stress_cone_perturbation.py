import math

import numpy as np
import pytest

from openai_ns_reconstruction.stress_cone import PositiveSignedStressDecomposition
from openai_ns_reconstruction.stress_cone_perturbation import (
    CovariancePerturbationCertificate,
)


def _reference() -> PositiveSignedStressDecomposition:
    return PositiveSignedStressDecomposition(
        a=2.75,
        b=1.6,
        scale_minus=0.8,
        scale_plus=1.25,
        m=3.2,
        t=0.45,
    )


def test_actual_matrix_bounds_and_positive_coefficients_crosscheck_independently():
    # Synthetic perturbations exercise only the finite-dimensional implication;
    # they are not claimed to be paper pulse columns.
    cert = CovariancePerturbationCertificate(
        reference=_reference(),
        error_minus=(0.03, -0.04),
        error_plus=(-0.02, 0.01),
        normalized_error_bound=math.nextafter(0.05, math.inf),
    )

    matrix = cert.actual_matrix
    target = cert.reference.target

    independent_det = float(np.linalg.det(matrix))
    independent_inverse_norm = float(np.linalg.norm(np.linalg.inv(matrix), ord=2))
    independent_coefficients = np.linalg.solve(matrix, target)
    reference_coefficients = cert.reference_squared_amplitudes

    assert cert.actual_determinant == pytest.approx(independent_det, rel=3e-14, abs=3e-14)
    assert cert.actual_determinant < 0.0
    assert abs(independent_det) >= cert.determinant_abs_lower_bound
    assert independent_inverse_norm <= cert.inverse_two_norm_upper

    assert cert.actual_squared_amplitudes == pytest.approx(
        independent_coefficients, rel=3e-14, abs=3e-14
    )
    independent_coefficient_error = float(
        np.linalg.norm(independent_coefficients - reference_coefficients, ord=2)
    )
    assert independent_coefficient_error <= cert.coefficient_error_two_norm_upper
    assert float(np.min(independent_coefficients)) >= cert.squared_amplitude_lower_bound
    assert np.all(cert.actual_amplitudes > 0.0)
    assert cert.reconstructed_target == pytest.approx(target, rel=3e-14, abs=3e-14)


def test_zero_perturbation_reduces_to_reference_matrix_and_solve():
    cert = CovariancePerturbationCertificate(
        reference=_reference(),
        error_minus=(0.0, 0.0),
        error_plus=(0.0, 0.0),
        normalized_error_bound=0.0,
    )

    assert cert.actual_matrix == pytest.approx(cert.reference.matrix)
    assert cert.actual_squared_amplitudes == pytest.approx(cert.reference_squared_amplitudes)
    assert cert.squared_amplitude_lower_bound > 0.0


def test_concrete_column_error_must_fit_the_certified_uniform_bound():
    with pytest.raises(ValueError, match="column error exceeds certified bound"):
        CovariancePerturbationCertificate(
            reference=_reference(),
            error_minus=(0.06, 0.0),
            error_plus=(0.0, 0.0),
            normalized_error_bound=0.05,
        )


def test_large_error_bound_fails_before_claiming_nondegeneracy():
    with pytest.raises(ValueError, match="determinant nondegeneracy"):
        CovariancePerturbationCertificate(
            reference=_reference(),
            error_minus=(0.0, 0.0),
            error_plus=(0.0, 0.0),
            normalized_error_bound=1.2,
        )


def test_intermediate_error_bound_can_keep_determinant_but_lose_positive_cone():
    with pytest.raises(ValueError, match="positive squared amplitudes"):
        CovariancePerturbationCertificate(
            reference=_reference(),
            error_minus=(0.0, 0.0),
            error_plus=(0.0, 0.0),
            normalized_error_bound=0.5,
        )


def test_invalid_error_data_fail_closed():
    with pytest.raises(ValueError, match=r"shape \(2,\)"):
        CovariancePerturbationCertificate(
            reference=_reference(),
            error_minus=(0.0,),
            error_plus=(0.0, 0.0),
            normalized_error_bound=0.1,
        )

    with pytest.raises(ValueError, match="finite"):
        CovariancePerturbationCertificate(
            reference=_reference(),
            error_minus=(math.nan, 0.0),
            error_plus=(0.0, 0.0),
            normalized_error_bound=0.1,
        )

    with pytest.raises(ValueError, match="nonnegative"):
        CovariancePerturbationCertificate(
            reference=_reference(),
            error_minus=(0.0, 0.0),
            error_plus=(0.0, 0.0),
            normalized_error_bound=-1e-6,
        )
