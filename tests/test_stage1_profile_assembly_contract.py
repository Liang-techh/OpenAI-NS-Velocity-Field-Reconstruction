from dataclasses import replace
from decimal import Decimal

import pytest

from openai_ns_reconstruction.axis_coefficient_reference_state import (
    WINDOW_LEFT,
    WINDOW_RIGHT,
    actual_schedule_reference_axis_state,
)
from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.stage1_profile_assembly_contract import (
    NATURAL_AMPLITUDE_CONVENTION,
    Stage1ProfileBackendIdentity,
    actual_schedule_profile_assembly_contract,
)
from openai_ns_reconstruction.stage1_scale_chain_wide import (
    diagnose_actual_schedule_scale_chain_wide,
)


@pytest.fixture(scope="module")
def schedule_data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


@pytest.fixture(scope="module")
def contract(schedule_data: TailData):
    return actual_schedule_profile_assembly_contract(schedule_data, 0.05)


def test_contract_is_owned_by_the_same_actual_schedule_reference_and_wide_scale(
    schedule_data: TailData,
    contract,
) -> None:
    reference = actual_schedule_reference_axis_state(schedule_data, 0.05)
    diagnostic = diagnose_actual_schedule_scale_chain_wide(schedule_data, 0.05)

    assert diagnostic.upstream_obstruction is None
    assert diagnostic.scale is not None
    assert contract.identity.schedule_key == (
        schedule_data.core.P,
        schedule_data.core.m,
        schedule_data.core.lam,
        schedule_data.core.wait,
        schedule_data.h,
    )
    assert contract.identity.j == 0.05
    assert contract.identity.sigma == reference.reference.sigma
    assert contract.identity.epsilon == reference.epsilon
    assert contract.identity.Lambda == diagnostic.scale.Lambda
    assert contract.identity.log_C == diagnostic.scale.C.exponent_upper
    assert contract.identity.window_left == WINDOW_LEFT
    assert contract.identity.window_right == WINDOW_RIGHT
    assert contract.identity.amplitude_convention == NATURAL_AMPLITUDE_CONVENTION
    assert contract.picard.Lambda == contract.identity.Lambda
    assert contract.paper_exact is False
    assert contract.full_reconstruction is False

    contract.require_backend_identity(contract.identity)


@pytest.mark.parametrize(
    ("field", "replacement"),
    (
        ("j", 0.049),
        ("sigma", 0.5),
        ("epsilon", 0.5),
        ("Lambda", Decimal(2)),
        ("log_C", Decimal(3)),
        ("window_left", -1.0),
        ("window_right", 1.0),
    ),
)
def test_contract_rejects_any_fixed_point_identity_drift(
    contract,
    field: str,
    replacement,
) -> None:
    candidate = replace(contract.identity, **{field: replacement})
    with pytest.raises(ValueError, match=field):
        contract.require_backend_identity(candidate)


def test_contract_rejects_a_different_actual_schedule(contract) -> None:
    changed = replace(
        contract.identity,
        schedule_key=(
            contract.identity.schedule_key[0] + 1.0,
            *contract.identity.schedule_key[1:],
        ),
    )
    with pytest.raises(ValueError, match="schedule_key"):
        contract.require_backend_identity(changed)


def test_identity_refuses_an_unpinned_amplitude_convention(contract) -> None:
    with pytest.raises(ValueError, match="amplitude convention"):
        replace(contract.identity, amplitude_convention="sampled amplitude")


def test_contract_requires_typed_backend_identity(contract) -> None:
    with pytest.raises(TypeError, match="Stage1ProfileBackendIdentity"):
        contract.require_backend_identity(object())  # type: ignore[arg-type]


def test_identity_requires_wide_decimal_scale(contract) -> None:
    with pytest.raises(ValueError, match="finite Decimal"):
        Stage1ProfileBackendIdentity(
            schedule_key=contract.identity.schedule_key,
            j=contract.identity.j,
            sigma=contract.identity.sigma,
            epsilon=contract.identity.epsilon,
            Lambda=1.0,  # type: ignore[arg-type]
            log_C=contract.identity.log_C,
            window_left=contract.identity.window_left,
            window_right=contract.identity.window_right,
        )