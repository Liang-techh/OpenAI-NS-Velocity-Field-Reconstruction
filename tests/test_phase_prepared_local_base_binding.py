from fractions import Fraction

import pytest

from openai_ns_reconstruction.phase_large_band_base_source import (
    LargeBandBaseSourceBinding,
    LargeBandBaseSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_local_base import (
    AsymptoticSlowLabel,
    LargeBandLocalBaseAdmission,
    LargeBandLocalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_physical_base import (
    LargeBandPhysicalBaseBinding,
    PrimaryGeometryPhysicalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_prepared_source import (
    LargeBandPreparedSourceBinding,
    PrimaryGeometryPreparedSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_scale import LargeBandPhaseScaleCertificate
from openai_ns_reconstruction.phase_prepared_local_base_binding import (
    PINNED_BASE_CHART_LOCAL_BASE_THEOREM,
    PINNED_LOCAL_BASE_BOUNDS,
    PINNED_PREPARED_BASE_FIELD,
    PreparedBaseLocalBoundsWitness,
    SlowBoxPreparedLocalBaseBinding,
)
from openai_ns_reconstruction.slow_base_field_binding import (
    SlowBoxBaseFieldBindingWitness,
    TheoremBackedSlowBoxBaseFieldBinding,
)


ELL = 100000
A = (7, -3, 11)
SOURCE_ID = "formal-base-provider"
SOURCE_REVISION = "fixture-r1"
PREPARED_ID = "prepared-fixture-001"
M = 2.0


def _physical_binding():
    scale = LargeBandPhaseScaleCertificate(ell=ELL, h=Fraction(1, 200), M=M)
    local = LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ELL, a=A, sigma=1),
        M=M,
        diameter_times_S3=0.75,
        second_F0_over_M=0.8,
        second_G0_over_M=0.7,
        first_F0_over_M=0.6,
        value_F_error_over_M_epsilon2=0.5,
        first_F_error_over_M_epsilon2=0.4,
        first_G_error_over_M_epsilon2=0.3,
        radial_F_over_M=0.9,
        axial_F_over_M=0.8,
        axial_G_over_M=0.7,
        evidence_kind="formal-theorem",
        provenance="synthetic theorem fixture for LocalBase admission",
        convex_domain_certified=True,
        differentiability_certified=True,
        representative_membership_certified=True,
        local_base_suprema_certified=True,
    )
    base_source = LargeBandBaseSourceBinding(
        admission=LargeBandLocalBaseAdmission(scale=scale, witness=local),
        source=LargeBandBaseSourceWitness(
            ell=ELL,
            a=A,
            M=M,
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for source identity",
            provider_contract_certified=True,
            sign_independent_source_certified=True,
            local_base_source_identity_certified=True,
        ),
    )
    prepared = LargeBandPreparedSourceBinding(
        base_source=base_source,
        prepared=PrimaryGeometryPreparedSourceWitness(
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            prepared_instance_id=PREPARED_ID,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for Prepared identity",
            exists_prepared_application_certified=True,
            prepared_choice_identity_certified=True,
            family_base_source_identity_certified=True,
            phases_reuse_prepared_certified=True,
            both_signs_share_prepared_certified=True,
        ),
    )
    return LargeBandPhysicalBaseBinding(
        prepared_source=prepared,
        witness=PrimaryGeometryPhysicalBaseWitness(
            ell=ELL,
            a=A,
            source_id=SOURCE_ID,
            source_revision=SOURCE_REVISION,
            prepared_instance_id=PREPARED_ID,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for physical-base chain",
            local_base_same_frequency_axial_certified=True,
            construction_frequency_certified=True,
            construction_axial_certified=True,
            carrier_time_positive_certified=True,
            prepared_radius_positive_certified=True,
            frequency_eq_physical_certified=True,
            axial_eq_physical_certified=True,
            both_signs_share_physical_base_certified=True,
        ),
    )


def _base_binding():
    physical = _physical_binding()
    interface = SlowBoxBaseFieldBindingWitness(
        box_key=physical.box_key,
        base_field_family_id="formal:background-family#section5",
        normalized_base_field_id="formal:normalized-base-field#beta",
        tangential_jet_provider_id="formal:tangential-jet-provider#beta",
        evaluation_application_id="formal:base-field-evaluation#beta",
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        prepared_instance_id=PREPARED_ID,
        producer_kind="formal-interface-export",
        provenance="synthetic interface metadata; not a materialized paper field",
        box_identity_certified=True,
        sign_independence_certified=True,
        normalized_chart_coordinates_certified=True,
        provider_derivation_certified=True,
        application_certified=True,
    )
    return TheoremBackedSlowBoxBaseFieldBinding(interface=interface, physical_base=physical)


def _witness(**changes):
    kwargs = dict(
        box_key=(ELL, A),
        source_id=SOURCE_ID,
        source_revision=SOURCE_REVISION,
        prepared_instance_id=PREPARED_ID,
        evidence_kind="formal-theorem",
        provenance="synthetic theorem application metadata for Prepared.base",
        prepared_base_application_certified=True,
        active_cell_index_identity_certified=True,
        frequency_axial_field_identity_certified=True,
        local_base_bounds_application_certified=True,
        epsilon_band_identity_certified=True,
        both_signs_share_local_base_certified=True,
    )
    kwargs.update(changes)
    return PreparedBaseLocalBoundsWitness(**kwargs)


def test_bridge_reuses_embedded_localbase_and_scale_without_truth_upgrade():
    base = _base_binding()
    linked = SlowBoxPreparedLocalBaseBinding(base=base, witness=_witness())

    admission = base.physical_base.prepared_source.base_source.admission
    assert linked.local_base_admission is admission
    assert linked.phase_scale is admission.scale
    assert linked.sign_pair == base.physical_base.prepared_source.sign_pair
    assert linked.simplified_phase_error_bound == admission.scale.simplified_phase_error_bound
    assert linked.simplified_rounded_normal_bound == admission.scale.simplified_rounded_normal_bound
    assert linked.prepared_base_local_bounds_linked is True
    assert linked.uniform_phase_scale_and_localbase_hypotheses_linked is True
    assert linked.status == "formal-structure"
    assert linked.actual_prepared_base_application_machine_replayed is False
    assert linked.actual_physical_base_values_materialized is False
    assert linked.uniform_eq_7_9_to_7_11_verified is False
    assert linked.paper_exact_velocity_available is False


@pytest.mark.parametrize(
    ("field", "value", "match"),
    [
        ("box_key", (ELL, (8, -3, 11)), "box key"),
        ("source_revision", "fixture-r2", "source identity"),
        ("prepared_instance_id", "different-prepared", "Prepared identity"),
    ],
)
def test_bridge_rejects_cross_wired_theorem_application(field, value, match):
    with pytest.raises(ValueError, match=match):
        SlowBoxPreparedLocalBaseBinding(
            base=_base_binding(),
            witness=_witness(**{field: value}),
        )


@pytest.mark.parametrize("evidence", ["sampled", "fitted", "numeric-scan", "paper-exact"])
def test_non_theorem_evidence_is_rejected(evidence):
    with pytest.raises(ValueError, match="analytic-theorem or formal-theorem"):
        _witness(evidence_kind=evidence)


@pytest.mark.parametrize(
    "field",
    [
        "prepared_base_application_certified",
        "active_cell_index_identity_certified",
        "frequency_axial_field_identity_certified",
        "local_base_bounds_application_certified",
        "epsilon_band_identity_certified",
        "both_signs_share_local_base_certified",
    ],
)
def test_missing_theorem_certificate_fails_closed(field):
    with pytest.raises(ValueError, match="theorem-certified true"):
        _witness(**{field: False})


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("prepared_base_field", "PrimaryGeometryAssembly.Prepared.fake"),
        ("local_base_bounds_decl", "PhaseEstimates.FakeBounds"),
        ("base_chart_local_base_theorem", "BaseChartJets.Estimates.fake"),
    ],
)
def test_formal_source_symbols_are_pinned(field, value):
    with pytest.raises(ValueError, match="pinned formal source"):
        _witness(**{field: value})


def test_pinned_symbols_match_the_formal_chain():
    witness = _witness()
    assert witness.prepared_base_field == PINNED_PREPARED_BASE_FIELD
    assert witness.local_base_bounds_decl == PINNED_LOCAL_BASE_BOUNDS
    assert witness.base_chart_local_base_theorem == PINNED_BASE_CHART_LOCAL_BASE_THEOREM


@pytest.mark.parametrize(
    "box_key",
    [
        (0, (1, 2, 3)),
        ("100000", (1, 2, 3)),
        (ELL, (1, 2)),
        (ELL, (1, True, 3)),
    ],
)
def test_malformed_box_keys_fail_closed(box_key):
    with pytest.raises(ValueError, match="box_key"):
        _witness(box_key=box_key)


def test_wrong_composition_types_are_rejected():
    with pytest.raises(TypeError, match="TheoremBackedSlowBoxBaseFieldBinding"):
        SlowBoxPreparedLocalBaseBinding(object(), _witness())
    with pytest.raises(TypeError, match="PreparedBaseLocalBoundsWitness"):
        SlowBoxPreparedLocalBaseBinding(_base_binding(), object())
