"""Fail-closed bridge from a bound slow box to ``Prepared.base`` LocalBase data.

The Section 6/7 slow-box binding already identifies one sign-free box with the
pinned ``PrimaryGeometryAssembly.Prepared`` object and with the corresponding
``FinalSlowBase.velocity`` components.  The pinned formal source also stores the
actual phase hypotheses in the field ``PrimaryGeometryAssembly.Prepared.base``:
for every active ``CellIndex`` it returns ``PhaseEstimates.LocalBaseBounds`` for
the same frequency/axial fields.

This module records that theorem application as a separate, fail-closed
certificate.  It deliberately reuses the ``LargeBandLocalBaseAdmission`` already
embedded in the physical-base binding, so no caller can inject a second set of
``M``, band, or derivative-error ratios.

No Lean object or field value is evaluated here.  A successful binding therefore
links theorem identities/hypotheses only; it does not machine-replay
``Prepared.base``, materialize the background, prove Eqs. (7.9)--(7.11), or
construct an oscillatory wave.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

from .phase_large_band_prepared_source import (
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_FILE,
    PINNED_LEAN_REPOSITORY,
)
from .slow_base_field_binding import TheoremBackedSlowBoxBaseFieldBinding


PINNED_PREPARED_BASE_FIELD = "PrimaryGeometryAssembly.Prepared.base"
PINNED_LOCAL_BASE_BOUNDS = "PhaseEstimates.LocalBaseBounds"
PINNED_BASE_CHART_LOCAL_BASE_THEOREM = "BaseChartJets.Estimates.localBaseBounds"

_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _box_key(value: object) -> tuple[int, tuple[int, int, int]]:
    if not isinstance(value, tuple) or len(value) != 2:
        raise ValueError("box_key must be (ell, a)")
    ell, a = value
    if isinstance(ell, bool) or not isinstance(ell, Integral) or int(ell) < 1:
        raise ValueError("box_key ell must be a positive integer")
    if not isinstance(a, tuple) or len(a) != 3:
        raise ValueError("box_key a must be a length-three integer tuple")
    if any(isinstance(x, bool) or not isinstance(x, Integral) for x in a):
        raise ValueError("box_key a must be a length-three integer tuple")
    return int(ell), tuple(int(x) for x in a)


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


def _pinned(value: object, expected: str, name: str) -> str:
    text = _nonempty_text(value, name)
    if text != expected:
        raise ValueError(f"{name} must match the pinned formal source exactly")
    return text


@dataclass(frozen=True)
class PreparedBaseLocalBoundsWitness:
    """Identity of one theorem application of ``Prepared.base``.

    The booleans are supplied by a theorem/export layer.  Python only checks
    that the application names the same box/source/Prepared object as the
    already-admitted physical-base chain and that the formal symbols are pinned.
    """

    box_key: tuple[int, tuple[int, int, int]]
    source_id: str
    source_revision: str
    prepared_instance_id: str
    evidence_kind: str
    provenance: str

    prepared_base_application_certified: bool
    active_cell_index_identity_certified: bool
    frequency_axial_field_identity_certified: bool
    local_base_bounds_application_certified: bool
    epsilon_band_identity_certified: bool
    both_signs_share_local_base_certified: bool

    lean_repository: str = PINNED_LEAN_REPOSITORY
    lean_commit: str = PINNED_LEAN_COMMIT
    lean_file: str = PINNED_LEAN_FILE
    prepared_base_field: str = PINNED_PREPARED_BASE_FIELD
    local_base_bounds_decl: str = PINNED_LOCAL_BASE_BOUNDS
    base_chart_local_base_theorem: str = PINNED_BASE_CHART_LOCAL_BASE_THEOREM

    def __post_init__(self) -> None:
        object.__setattr__(self, "box_key", _box_key(self.box_key))
        for name in ("source_id", "source_revision", "prepared_instance_id", "provenance"):
            object.__setattr__(self, name, _nonempty_text(getattr(self, name), name))

        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; "
                "sampled/fitted/numeric-scan evidence is rejected"
            )

        for name in (
            "prepared_base_application_certified",
            "active_cell_index_identity_certified",
            "frequency_axial_field_identity_certified",
            "local_base_bounds_application_certified",
            "epsilon_band_identity_certified",
            "both_signs_share_local_base_certified",
        ):
            _strict_true(getattr(self, name), name)

        pins = (
            ("lean_repository", PINNED_LEAN_REPOSITORY),
            ("lean_commit", PINNED_LEAN_COMMIT),
            ("lean_file", PINNED_LEAN_FILE),
            ("prepared_base_field", PINNED_PREPARED_BASE_FIELD),
            ("local_base_bounds_decl", PINNED_LOCAL_BASE_BOUNDS),
            ("base_chart_local_base_theorem", PINNED_BASE_CHART_LOCAL_BASE_THEOREM),
        )
        for name, expected in pins:
            object.__setattr__(self, name, _pinned(getattr(self, name), expected, name))

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def prepared_key(self) -> tuple[str, str, str]:
        return self.source_id, self.source_revision, self.prepared_instance_id


@dataclass(frozen=True)
class SlowBoxPreparedLocalBaseBinding:
    """Compose the provider/physical-base identity with the pinned LocalBase field.

    The quantitative LocalBase witness and exact large-band scale are not passed
    independently: they are read from
    ``base.physical_base.prepared_source.base_source.admission``.  This prevents
    a theorem application for one box from being combined with LocalBase ratios
    or scale data from another box.
    """

    base: TheoremBackedSlowBoxBaseFieldBinding
    witness: PreparedBaseLocalBoundsWitness

    def __post_init__(self) -> None:
        if not isinstance(self.base, TheoremBackedSlowBoxBaseFieldBinding):
            raise TypeError("base must be a TheoremBackedSlowBoxBaseFieldBinding")
        if not isinstance(self.witness, PreparedBaseLocalBoundsWitness):
            raise TypeError("witness must be a PreparedBaseLocalBoundsWitness")

        admission = self.base.physical_base.prepared_source.base_source.admission
        if self.witness.box_key != self.base.box_key:
            raise ValueError("Prepared.base box key must match the bound slow box")
        if self.witness.source_key != self.base.source_key:
            raise ValueError("Prepared.base source identity must match the bound base source")
        if self.witness.prepared_key != self.base.prepared_key:
            raise ValueError("Prepared.base Prepared identity must match the selected Prepared source")
        if admission.witness.label.box_key != self.base.box_key:
            raise ValueError("embedded LocalBase admission must use the bound slow box")
        if not all(admission.admission_checks().values()):
            raise ValueError("embedded LocalBase admission must remain fully certified")
        if not self.base.theorem_physical_base_identity_linked:
            raise ValueError("physical-base theorem identity link must remain certified")

    @property
    def local_base_admission(self):
        return self.base.physical_base.prepared_source.base_source.admission

    @property
    def phase_scale(self):
        return self.local_base_admission.scale

    @property
    def sign_pair(self):
        return self.base.physical_base.prepared_source.sign_pair

    @property
    def simplified_phase_error_bound(self) -> float:
        return self.phase_scale.simplified_phase_error_bound

    @property
    def simplified_rounded_normal_bound(self) -> float:
        return self.phase_scale.simplified_rounded_normal_bound

    @property
    def prepared_base_local_bounds_linked(self) -> bool:
        return True

    @property
    def uniform_phase_scale_and_localbase_hypotheses_linked(self) -> bool:
        return True

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def actual_prepared_base_application_machine_replayed(self) -> bool:
        return False

    @property
    def actual_physical_base_values_materialized(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
