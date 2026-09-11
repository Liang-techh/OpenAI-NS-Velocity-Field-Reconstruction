"""Wide Stage-1 scale bridge for the actual SchedulePressure construction.

This module starts from the existing current-main schedule/analytic-input chain,
reuses its componentwise chi-resolvent certificate, and changes only the two
scalar layers that currently exceed binary64: controlled-remainder propagation
and the final Lambda/C selection.

It does not solve the coefficient-space fixed point.  Reaching a symbolic C
threshold is therefore a representation-layer advance inside the existing
formal-structure boundary, not a paper-exact leading profile.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .axis_remainder_wide_bounds import (
    NaturalRemainderWideCertificate,
    natural_remainder_bound_certificate_wide,
)
from .natural_scale_selection_wide import (
    WideNaturalScaleSelection,
    select_natural_scale_wide,
)
from .outgoing_tail import TailData
from .stage1_scale_chain import (
    ActualScheduleScaleChainDiagnostic,
    diagnose_actual_schedule_scale_chain,
)


@dataclass(frozen=True)
class ActualScheduleWideScaleChainDiagnostic:
    """Actual-schedule path continued through wide remainder and Lambda/C."""

    narrow_prefix: ActualScheduleScaleChainDiagnostic
    remainder: Optional[NaturalRemainderWideCertificate]
    scale: Optional[WideNaturalScaleSelection]
    upstream_obstruction: Optional[str]

    @property
    def reached_wide_remainder(self) -> bool:
        return self.remainder is not None

    @property
    def reached_wide_scale_selection(self) -> bool:
        return self.scale is not None


def diagnose_actual_schedule_scale_chain_wide(
    data: TailData,
    j: float,
) -> ActualScheduleWideScaleChainDiagnostic:
    """Continue the landed actual-schedule chain without narrowing wide bounds.

    The narrow prefix remains authoritative for the schedule, common analytic
    neighborhood, per-field Cauchy bounds, and natural-resolvent majorant.  If
    that prefix has not reached a representable ``analytic_norms`` object, this
    adapter stops rather than fabricating a substitute.  Otherwise it evaluates
    exactly the same controlled-remainder ledger in wide Decimal arithmetic and
    carries those wide bounds into the pinned Lambda/C selection.
    """

    prefix = diagnose_actual_schedule_scale_chain(data, j)
    analytic = prefix.analytic_norms
    if analytic is None:
        detail = "actual-schedule prefix did not reach representable analytic operator norms"
        if prefix.obstruction is not None:
            detail = (
                "actual-schedule factorial majorant remains outside the narrow prefix: "
                f"term {prefix.obstruction.term_index}"
            )
        elif prefix.propagation_obstruction is not None:
            detail = (
                f"actual-schedule prefix stopped at "
                f"{prefix.propagation_obstruction.stage}: "
                f"{prefix.propagation_obstruction.detail}"
            )
        return ActualScheduleWideScaleChainDiagnostic(
            narrow_prefix=prefix,
            remainder=None,
            scale=None,
            upstream_obstruction=detail,
        )

    try:
        operators = analytic.operator_norm_bounds()
    except Exception as exc:
        return ActualScheduleWideScaleChainDiagnostic(
            narrow_prefix=prefix,
            remainder=None,
            scale=None,
            upstream_obstruction=f"operator norm conversion failed: {exc}",
        )

    axis_data = analytic.axis_data_norm_bounds(h=data.h)
    remainder = natural_remainder_bound_certificate_wide(
        operators=operators,
        data=axis_data,
        amplitude_norm_upper=analytic.amplitude_norm_upper,
    )
    scale = select_natural_scale_wide(
        remainder_bound=remainder.remainder_bound_upper,
        remainder_lipschitz=remainder.remainder_lipschitz_upper,
        phase_real_part_sup=(
            prefix.schedule_inputs.neighborhood.axis_phase_real_part_sup_upper
        ),
    )

    return ActualScheduleWideScaleChainDiagnostic(
        narrow_prefix=prefix,
        remainder=remainder,
        scale=scale,
        upstream_obstruction=None,
    )
