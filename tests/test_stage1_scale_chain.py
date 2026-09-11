from decimal import Decimal
import math

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


def test_actual_schedule_chain_uses_chi_specific_resolvent_bound() -> None:
    result = diagnose_actual_schedule_scale_chain(_data(), 0.05)

    assert result.schedule_inputs.neighborhood.radius > 0.0
    assert result.coefficient_family_norm_upper > 0.0
    assert result.resolvent_majorant_parameter_upper > 0.0

    # The old common-B route let the much larger gradient bound contaminate chi
    # and produced a positive factorial-majorant term above binary64.  The
    # theorem-faithful componentwise boundedAxisElement ledger removes that
    # artificial coupling: K uses only 12 * B_chi.
    assert result.obstruction is None
    assert result.analytic_norms is not None
    assert result.resolvent_majorant_parameter_upper < 100_000.0
    assert (
        result.resolvent_majorant_parameter_upper
        < 2560.0 * result.coefficient_family_norm_upper
    )
    assert math.isfinite(result.analytic_norms.resolvent_majorant.series_upper)

    # Fail closed at the next actual representability boundary.  This is a
    # downstream conservative-bound overflow, not evidence that the theorem's
    # true fixed point or resolvent fails to exist.
    assert result.propagation_obstruction is not None
    assert result.propagation_obstruction.stage == "remainder-propagation"
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
    # This standalone large-K regression keeps coverage of the Decimal witness
    # even though the actual schedule no longer needs the common-B K bound.
    witness = first_binary64_majorant_obstruction(1.0e14, max_terms=64)
    assert witness is not None
    assert isinstance(witness.term_lower, Decimal)
    assert witness.term_lower > witness.binary64_max
