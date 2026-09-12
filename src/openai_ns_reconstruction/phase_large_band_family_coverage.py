"""Fail-closed finite-scope coverage for Section 7 large-band theorem witnesses.

The per-family gates in :mod:`phase_large_band_phase_estimates`,
:mod:`phase_large_band_coordinate_errors`, and :mod:`phase_large_band_damping`
are deliberately local: each proves only that one signed family label is tied
to the pinned theorem metadata and the same sign-free base-source revision.
Those local gates are not, by themselves, evidence that *every* active family
in a claimed scope was covered.

This module adds only that missing bookkeeping layer.  A caller must provide a
theorem-provenanced declaration of the exact finite active signed family set
for one large band and one base-source revision.  Coverage is admitted only if
phase, coordinate, and damping admissions each cover exactly that declared
set, without duplicates or extras, and all three paths reuse the same upstream
family admission for every label.

The declared scope itself is still caller supplied: Python does not derive the
paper's active carrier/slot set or verify the theorem proving that the tuple is
exact.  Therefore successful coverage is ``formal-structure`` relative to the
scope witness.  It does not establish global all-band uniformity, Eqs.
(7.9)--(7.11), or paper-exact velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .phase_large_band_coordinate_errors import LargeBandCoordinateErrorsAdmission
from .phase_large_band_damping import LargeBandDampingErrorAdmission
from .phase_large_band_local_base import AsymptoticSlowLabel
from .phase_large_band_phase_estimates import LargeBandPhaseEstimatesAdmission


_ACCEPTED_EVIDENCE = frozenset({"analytic-theorem", "formal-theorem"})


def _nonempty_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value.strip()


def _strict_true(value: object, name: str) -> bool:
    if value is not True:
        raise ValueError(f"{name} must be theorem-certified true")
    return True


def _unique_label_map(items: tuple[object, ...], attr: str, name: str) -> dict[AsymptoticSlowLabel, object]:
    result: dict[AsymptoticSlowLabel, object] = {}
    for item in items:
        label = getattr(item, attr)
        if label in result:
            raise ValueError(f"duplicate {name} admission for signed label {label!r}")
        result[label] = item
    return result


@dataclass(frozen=True)
class LargeBandActiveFamilyScopeWitness:
    """Theorem-facing declaration of one exact finite active signed scope.

    ``labels`` is not discovered here.  The certification flags must come from
    analytic/formal theorem provenance establishing that, for ``scope_id``, the
    tuple is the complete active carrier/slot family at this band and source
    revision.  Mechanical sign-pair closure is checked in addition to that
    provenance claim.
    """

    scope_id: str
    ell: int
    source_id: str
    source_revision: str
    labels: tuple[AsymptoticSlowLabel, ...]
    evidence_kind: str
    provenance: str
    finite_scope_certified: bool
    exact_active_signed_family_certified: bool
    carrier_slot_scope_certified: bool
    sign_duplication_certified: bool

    def __post_init__(self) -> None:
        object.__setattr__(self, "scope_id", _nonempty_text(self.scope_id, "scope_id"))
        if isinstance(self.ell, bool) or not isinstance(self.ell, int) or self.ell <= 0:
            raise ValueError("ell must be a positive integer")
        object.__setattr__(self, "source_id", _nonempty_text(self.source_id, "source_id"))
        object.__setattr__(
            self, "source_revision", _nonempty_text(self.source_revision, "source_revision")
        )
        object.__setattr__(self, "provenance", _nonempty_text(self.provenance, "scope provenance"))
        if self.evidence_kind not in _ACCEPTED_EVIDENCE:
            raise ValueError(
                "evidence_kind must be analytic-theorem or formal-theorem; sampled/fitted evidence is rejected"
            )
        if not isinstance(self.labels, tuple) or not self.labels:
            raise ValueError("labels must be a nonempty tuple of active signed labels")
        if any(not isinstance(label, AsymptoticSlowLabel) for label in self.labels):
            raise TypeError("every scope label must be an AsymptoticSlowLabel")
        if any(label.ell != self.ell for label in self.labels):
            raise ValueError("every scope label must use the declared ell")
        if len(set(self.labels)) != len(self.labels):
            raise ValueError("scope labels must be unique")

        labels = set(self.labels)
        for label in self.labels:
            if label.opposite() not in labels:
                raise ValueError("scope labels must be closed under the Section 6 sign pair")

        for name in (
            "finite_scope_certified",
            "exact_active_signed_family_certified",
            "carrier_slot_scope_certified",
            "sign_duplication_certified",
        ):
            _strict_true(getattr(self, name), name)

    @property
    def expected_labels(self) -> frozenset[AsymptoticSlowLabel]:
        return frozenset(self.labels)

    @property
    def source_key(self) -> tuple[str, str]:
        return self.source_id, self.source_revision

    @property
    def box_keys(self) -> frozenset[tuple[int, tuple[int, int, int]]]:
        return frozenset(label.box_key for label in self.labels)


@dataclass(frozen=True)
class LargeBandFamilyEstimateCoverage:
    """Exact-set coverage of the admitted Section 7 theorem paths.

    The three admission tuples must each equal ``scope.expected_labels`` as a
    set.  For every signed label, the phase path used by coordinate errors and
    the family path used by damping must be structurally identical to the phase
    admission supplied for that label.  This prevents sparse cherry-picking or
    mixing source revisions while reporting a finite-scope uniform envelope.
    """

    scope: LargeBandActiveFamilyScopeWitness
    phase_admissions: tuple[LargeBandPhaseEstimatesAdmission, ...]
    coordinate_admissions: tuple[LargeBandCoordinateErrorsAdmission, ...]
    damping_admissions: tuple[LargeBandDampingErrorAdmission, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.scope, LargeBandActiveFamilyScopeWitness):
            raise TypeError("scope must be a LargeBandActiveFamilyScopeWitness")
        for seq, cls, name in (
            (self.phase_admissions, LargeBandPhaseEstimatesAdmission, "phase"),
            (self.coordinate_admissions, LargeBandCoordinateErrorsAdmission, "coordinate"),
            (self.damping_admissions, LargeBandDampingErrorAdmission, "damping"),
        ):
            if not isinstance(seq, tuple):
                raise TypeError(f"{name}_admissions must be a tuple")
            if any(not isinstance(item, cls) for item in seq):
                raise TypeError(f"every {name} admission must have type {cls.__name__}")

        phase = _unique_label_map(self.phase_admissions, "label", "phase")
        coordinate = _unique_label_map(self.coordinate_admissions, "label", "coordinate")
        damping = _unique_label_map(self.damping_admissions, "label", "damping")
        expected = set(self.scope.expected_labels)
        for name, mapping in (("phase", phase), ("coordinate", coordinate), ("damping", damping)):
            actual = set(mapping)
            if actual != expected:
                missing = expected - actual
                extra = actual - expected
                raise ValueError(
                    f"{name} coverage must equal the declared active scope; "
                    f"missing={len(missing)}, extra={len(extra)}"
                )

        for label in self.scope.labels:
            p = phase[label]
            c = coordinate[label]
            d = damping[label]
            if p.source_key != self.scope.source_key:
                raise ValueError("phase source revision must match the declared active scope")
            if c.source_key != self.scope.source_key:
                raise ValueError("coordinate source revision must match the declared active scope")
            if d.source_key != self.scope.source_key:
                raise ValueError("damping source revision must match the declared active scope")
            if p.family != c.phase.family:
                raise ValueError("coordinate path must reuse the same upstream family admission")
            if p.family != d.family:
                raise ValueError("damping path must reuse the same upstream family admission")
            if not all(p.admission_checks().values()):
                raise ValueError("phase admission is no longer internally certified")
            if not all(c.admission_checks().values()):
                raise ValueError("coordinate admission is no longer internally certified")
            if not all(d.admission_checks().values()):
                raise ValueError("damping admission is no longer internally certified")

        bounds = self.uniform_scope_bounds()
        if any(not math.isfinite(value) or value <= 0.0 for value in bounds.values()):
            raise ValueError("finite-scope theorem envelopes must be finite and positive")

    def uniform_scope_bounds(self) -> dict[str, float]:
        """Return maxima over the exact admitted finite scope.

        These are theorem-envelope maxima, not measured field errors.  They are
        uniform only over the supplied theorem-certified finite scope.
        """

        return {
            "phase_normal_error": max(item.normal_error_bound for item in self.phase_admissions),
            "phase_velocity_error": max(item.phase_velocity_bound for item in self.phase_admissions),
            "coordinate_error": max(
                max(
                    item.coordinate_error_A_bound,
                    item.coordinate_error_B_bound,
                    item.coordinate_error_C_bound,
                )
                for item in self.coordinate_admissions
            ),
            "damping_error": max(item.damping_theorem_envelope for item in self.damping_admissions),
        }

    @property
    def covered_labels(self) -> frozenset[AsymptoticSlowLabel]:
        return frozenset(item.label for item in self.phase_admissions)

    @property
    def scope_coverage_complete(self) -> bool:
        return self.covered_labels == self.scope.expected_labels

    @property
    def status(self) -> str:
        return "formal-structure"

    @property
    def exact_active_family_scope_machine_verified(self) -> bool:
        return False

    @property
    def theorem_witnesses_machine_verified(self) -> bool:
        return False

    @property
    def global_all_band_uniformity_verified(self) -> bool:
        return False

    @property
    def uniform_eq_7_9_to_7_11_verified(self) -> bool:
        return False

    @property
    def paper_exact_velocity_available(self) -> bool:
        return False
