import pytest

from openai_ns_reconstruction.phase_large_band_family_coverage import (
    LargeBandActiveFamilyScopeWitness,
)
from openai_ns_reconstruction.phase_large_band_local_base import AsymptoticSlowLabel
from openai_ns_reconstruction.phase_large_band_phase_estimates import (
    PINNED_LEAN_COMMIT,
    PINNED_LEAN_REPOSITORY,
)
from openai_ns_reconstruction.phase_large_band_primary_geometry_scope import (
    PINNED_CELL_DOMAIN_DEFINITION,
    PINNED_CELL_INDEX_DEFINITION,
    PINNED_FAMILY_DEFINITION,
    LargeBandPrimaryGeometryScopeAdmission,
    PrimaryGeometryBandSliceWitness,
)


ELL = 100000
A0 = (7, -3, 11)
A1 = (8, -3, 11)
SOURCE_ID = "formal-base-provider"
SOURCE_REVISION = "fixture-r1"


def _labels(boxes=(A0, A1)):
    return tuple(
        AsymptoticSlowLabel(ell=ELL, a=box, sigma=sigma)
        for box in boxes
        for sigma in (-1, 1)
    )


def _scope(**overrides):
    kwargs = dict(
        scope_id="fixture-primary-geometry-band-slice",
        ell=ELL,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        labels=_labels(),
        evidence_kind="formal-theorem",
        provenance="synthetic formal scope fixture; not manuscript data",
        finite_scope_certified=True,
        exact_active_signed_family_certified=True,
        carrier_slot_scope_certified=True,
        sign_duplication_certified=True,
    )
    kwargs.update(overrides)
    return LargeBandActiveFamilyScopeWitness(**kwargs)


def _witness(**overrides):
    kwargs = dict(
        scope_id="fixture-primary-geometry-band-slice",
        ell=ELL,
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        unsigned_boxes=(A0, A1),
        lean_repository=PINNED_LEAN_REPOSITORY,
        lean_commit=PINNED_LEAN_COMMIT,
        cell_index_definition=PINNED_CELL_INDEX_DEFINITION,
        cell_domain_definition=PINNED_CELL_DOMAIN_DEFINITION,
        family_definition=PINNED_FAMILY_DEFINITION,
        evidence_kind="formal-theorem",
        provenance="synthetic theorem export fixture; not a CellIndex enumeration",
        complete_unsigned_box_manifest_certified=True,
        cell_index_band_slice_certified=True,
        band_and_cell_projection_certified=True,
        cell_domain_carrier_identity_certified=True,
        family_index_identity_certified=True,
        fin2_sign_duplication_certified=True,
    )
    kwargs.update(overrides)
    return PrimaryGeometryBandSliceWitness(**kwargs)


def test_exact_unsigned_manifest_reconstructs_declared_signed_scope():
    admission = LargeBandPrimaryGeometryScopeAdmission(scope=_scope(), witness=_witness())

    assert admission.manifest_matches_scope is True
    assert admission.unsigned_box_count == 2
    assert admission.signed_family_count == 4
    assert admission.witness.expected_signed_labels == frozenset(_labels())


def test_missing_or_extra_unsigned_box_fails_exact_set_check():
    with pytest.raises(ValueError, match="unsigned box set must equal"):
        LargeBandPrimaryGeometryScopeAdmission(
            scope=_scope(),
            witness=_witness(unsigned_boxes=(A0,)),
        )

    extra = (9, -3, 11)
    with pytest.raises(ValueError, match="unsigned box set must equal"):
        LargeBandPrimaryGeometryScopeAdmission(
            scope=_scope(),
            witness=_witness(unsigned_boxes=(A0, A1, extra)),
        )


def test_scope_identity_and_base_source_revision_cannot_drift():
    with pytest.raises(ValueError, match="scope_id"):
        LargeBandPrimaryGeometryScopeAdmission(
            scope=_scope(scope_id="other-scope"),
            witness=_witness(),
        )

    with pytest.raises(ValueError, match="base source revision"):
        LargeBandPrimaryGeometryScopeAdmission(
            scope=_scope(source_revision="different-base-revision"),
            witness=_witness(),
        )


def test_formal_source_pins_fail_closed():
    for field, value, message in (
        ("lean_repository", "other/repo", "lean_repository"),
        ("lean_commit", "0" * 40, "lean_commit"),
        ("cell_index_definition", "Fake.CellIndex", "cell_index_definition"),
        ("cell_domain_definition", "Fake.cellDomain", "cell_domain_definition"),
        ("family_definition", "Fake.family", "family_definition"),
    ):
        with pytest.raises(ValueError, match=message):
            _witness(**{field: value})


def test_manifest_requires_theorem_evidence_and_all_identity_certificates():
    for kind in ("sampled", "fitted", "numeric-scan", ""):
        with pytest.raises(ValueError, match="evidence_kind"):
            _witness(evidence_kind=kind)

    for field in (
        "complete_unsigned_box_manifest_certified",
        "cell_index_band_slice_certified",
        "band_and_cell_projection_certified",
        "cell_domain_carrier_identity_certified",
        "family_index_identity_certified",
        "fin2_sign_duplication_certified",
    ):
        with pytest.raises(ValueError, match=field):
            _witness(**{field: False})


def test_duplicate_or_malformed_unsigned_boxes_are_rejected():
    with pytest.raises(ValueError, match="unsigned_boxes must be unique"):
        _witness(unsigned_boxes=(A0, A0))

    with pytest.raises(ValueError, match="length-three integer tuple"):
        _witness(unsigned_boxes=((1, 2),))


def test_success_keeps_actual_enumeration_and_paper_truth_false():
    admission = LargeBandPrimaryGeometryScopeAdmission(scope=_scope(), witness=_witness())

    assert admission.status == "formal-structure"
    assert admission.actual_cell_index_enumeration_machine_verified is False
    assert admission.actual_active_family_scope_machine_verified is False
    assert admission.actual_base_fields_verified is False
    assert admission.uniform_eq_7_9_to_7_11_verified is False
    assert admission.paper_exact_velocity_available is False
