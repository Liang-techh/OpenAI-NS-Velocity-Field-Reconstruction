from decimal import Decimal

import pytest

from openai_ns_reconstruction.outgoing_tail import OutgoingCoreParameters, TailData
from openai_ns_reconstruction.stage1_scale_chain import (
    diagnose_actual_schedule_scale_chain,
    first_binary64_majorant_obstruction,
)


def _data() -> TailData:
    return TailData(
        OutgoingCoreParameters(P=2.0, m=1.0, lam=0.05, wait=30.0),
        h=0.01,
    )


def test_actual_schedule_chain_identifies_first_binary64_certificate_obstruction() -> None:
    result = diagnose_actual_schedule_scale_chain(_data(), 0.05)

    assert result.schedule_inputs.neighborhood.radius > 0.0
    assert result.coefficient_family_norm_upper > 0.0
    assert result.resolvent_majorant_parameter_upper > 0.0

    witness = result.obstruction
    assert witness is not None
    assert witness.term_index <= 64
    assert witness.term_lower > witness.binary64_max
    assert witness.decimal_order_lower >= 308

    # Fail closed: once the conservative resolvent majorant does not fit the
    # current numeric representation, no remainder/Lambda/C value is fabricated.
    assert result.analytic_norms is None
    assert result.remainder is None
    assert result.scale is None
    assert result.reached_remainder_chain is False
    assert result.reached_scale_selection is False


def test_obstruction_is_about_the_positive_majorant_not_an_actual_resolvent_lower_bound() -> None:
    # For a modest K no finite-prefix overflow witness exists.  None is
    # intentionally inconclusive rather than being promoted to a fit certificate.
    assert first_binary64_majorant_obstruction(1.0, max_terms=64) is None

    with pytest.raises(ValueError, match="positive integer"):
        first_binary64_majorant_obstruction(1.0, max_terms=0)
    with pytest.raises(ValueError, match="finite and positive"):
        first_binary64_majorant_obstruction(0.0)


def test_witness_decimal_value_is_strictly_larger_than_binary64_max() -> None:
    # This standalone large-K regression checks the Decimal comparison itself;
    # the actual-schedule test above derives K from the theorem-side schedule.
    witness = first_binary64_majorant_obstruction(1.0e14, max_terms=64)
    assert witness is not None
    assert isinstance(witness.term_lower, Decimal)
    assert witness.term_lower > witness.binary64_max
