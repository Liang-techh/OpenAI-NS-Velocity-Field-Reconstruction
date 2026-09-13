from fractions import Fraction
import math

import pytest

from openai_ns_reconstruction.charts import DyadicChart
from openai_ns_reconstruction.phase import TangentialBaseJet
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
    PINNED_FINAL_SLOW_BASE_VELOCITY,
    LargeBandPhysicalBaseBinding,
    PrimaryGeometryPhysicalBaseWitness,
)
from openai_ns_reconstruction.phase_large_band_prepared_source import (
    LargeBandPreparedSourceBinding,
    PrimaryGeometryPreparedSourceWitness,
)
from openai_ns_reconstruction.phase_large_band_scale import (
    LargeBandPhaseScaleCertificate,
)
from openai_ns_reconstruction.slow_base_field_binding import (
    BoundFrozenLabelPhaseData,
    SlowBoxBaseFieldBindingWitness,
    TheoremBackedSlowBoxBaseFieldBinding,
    freeze_bound_label_phase_data,
)
from openai_ns_reconstruction.slow_labels import (
    ActiveSlowRepresentative,
    SlowBoxEnclosure,
    SlowLabel,
)


def _active_representative(*, sigma=1):
    chart = DyadicChart(ell=200, h=0.005)
    ratio = 1.25
    R0, Z0 = 1.0, 0.2
    T0 = ratio - Z0**2 * ratio ** (2.0 * chart.h)
    label = SlowLabel(chart.ell, (7, -3, 11), sigma)
    mesh = chart.S_star ** -3
    enclosure = SlowBoxEnclosure(
        chart=chart,
        label=label,
        center=(R0, Z0, T0),
        half_width=(mesh / 2.0,) * 3,
    )
    return ActiveSlowRepresentative(
        chart=chart,
        label=label,
        R0=R0,
        Z0=Z0,
        T0=T0,
        q_big=1.5 * chart.Q,
        X_a=0.35,
        X_b=0.45,
        enclosure=enclosure,
    )


def _binding(box_key=None, **changes):
    rep = _active_representative()
    kwargs = dict(
        box_key=rep.label.box_key if box_key is None else box_key,
        base_field_family_id="formal:background-family#section5",
        normalized_base_field_id="formal:normalized-base-field#beta",
        tangential_jet_provider_id="formal:tangential-jet-provider#beta",
        evaluation_application_id="formal:base-field-evaluation#beta",
        source_id="formal-base-provider",
        source_revision="fixture-r1",
        prepared_instance_id="prepared-fixture-001",
        producer_kind="formal-interface-export",
        provenance="synthetic interface metadata; not a materialized paper field",
        box_identity_certified=True,
        sign_independence_certified=True,
        normalized_chart_coordinates_certified=True,
        provider_derivation_certified=True,
        application_certified=True,
    )
    kwargs.update(changes)
    return SlowBoxBaseFieldBindingWitness(**kwargs)


class Provider:
    base_field_provider_id = "formal:tangential-jet-provider#beta"

    def __init__(self):
        self.calls = []

    def tangential_jet(self, *, chart, R, Z, T):
        self.calls.append((chart.ell, R, Z, T))
        return TangentialBaseJet(
            F=1.0,
            G=0.25,
            F_R=-3.0,
            F_Z=0.1,
            F_T=-0.2,
            G_R=0.0,
            G_Z=-0.1,
            G_T=0.3,
        )


def _physical_binding():
    ell = 100000
    a = (7, -3, 11)
    M = 2.0
    source_id = "formal-base-provider"
    source_revision = "fixture-r1"
    prepared_id = "prepared-fixture-001"

    scale = LargeBandPhaseScaleCertificate(
        ell=ell, h=Fraction(1, 200), M=M
    )
    local = LargeBandLocalBaseWitness(
        label=AsymptoticSlowLabel(ell=ell, a=a, sigma=1),
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
    base = LargeBandBaseSourceBinding(
        admission=LargeBandLocalBaseAdmission(scale=scale, witness=local),
        source=LargeBandBaseSourceWitness(
            ell=ell,
            a=a,
            M=M,
            source_id=source_id,
            source_revision=source_revision,
            evidence_kind="formal-theorem",
            provenance="synthetic theorem fixture for source identity",
            provider_contract_certified=True,
            sign_independent_source_certified=True,
            local_base_source_identity_certified=True,
        ),
    )
    prepared = LargeBandPreparedSourceBinding(
        base_source=base,
        prepared=PrimaryGeometryPreparedSourceWitness(
            source_id=source_id,
            source_revision=source_revision,
            prepared_instance_id=prepared_id,
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
            ell=ell,
            a=a,
            source_id=source_id,
            source_revision=source_revision,
            prepared_instance_id=prepared_id,
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


def test_bound_freeze_connects_sign_free_box_to_identified_provider_without_truth_upgrade():
    rep = _active_representative(sigma=1)
    provider = Provider()
    binding = _binding(rep.label.box_key)

    bound = freeze_bound_label_phase_data(
        rep,
        provider,
        binding,
        rectangle_radius=0.01,
        u_star=0.25,
    )

    assert isinstance(bound, BoundFrozenLabelPhaseData)
    assert provider.calls == [(rep.chart.ell, rep.R0, rep.Z0, rep.T0)]
    assert bound.frozen.source is rep
    assert bound.frozen.frame.lambda0 == pytest.approx(math.sqrt(2.0))
    assert bound.binding.box_key == rep.label.box_key
    assert bound.slow_box_bound_to_base_field_interface is True
    assert bound.label_signs_share_box_binding is True
    assert bound.status == "formal-structure"
    assert bound.base_field_construction_machine_replayed is False
    assert bound.paper_exact_base_fields_materialized is False
    assert bound.base_jet_values_machine_certified is False
    assert bound.uniform_phase_estimates_proved is False
    assert bound.paper_exact_velocity_available is False


def test_opposite_sign_uses_the_same_sign_free_base_field_binding():
    rep = _active_representative(sigma=1)
    opposite_label = rep.label.opposite()
    opposite = ActiveSlowRepresentative(
        chart=rep.chart,
        label=opposite_label,
        R0=rep.R0,
        Z0=rep.Z0,
        T0=rep.T0,
        q_big=rep.q_big,
        X_a=rep.X_a,
        X_b=rep.X_b,
        enclosure=rep.enclosure,
    )
    binding = _binding(rep.label.box_key)

    plus = freeze_bound_label_phase_data(
        rep, Provider(), binding, rectangle_radius=0.01, u_star=0.25
    )
    minus = freeze_bound_label_phase_data(
        opposite, Provider(), binding, rectangle_radius=0.01, u_star=0.25
    )

    assert plus.binding == minus.binding
    assert plus.frozen.source.label.sigma == 1
    assert minus.frozen.source.label.sigma == -1
    assert plus.frozen.source.label.box_key == minus.frozen.source.label.box_key


def test_cross_wired_box_is_rejected_before_provider_evaluation():
    rep = _active_representative()
    provider = Provider()
    wrong = _binding((rep.chart.ell, (8, -3, 11)))

    with pytest.raises(ValueError, match="identical bound slow box"):
        freeze_bound_label_phase_data(
            rep, provider, wrong, rectangle_radius=0.01, u_star=0.25
        )
    assert provider.calls == []


def test_cross_wired_provider_identity_is_rejected_before_evaluation():
    rep = _active_representative()
    provider = Provider()
    provider.base_field_provider_id = "formal:other-provider"
    with pytest.raises(ValueError, match="base_field_provider_id"):
        freeze_bound_label_phase_data(
            rep, provider, _binding(), rectangle_radius=0.01, u_star=0.25
        )
    assert provider.calls == []


def test_theorem_bridge_connects_provider_identity_to_pinned_final_slow_base():
    physical = _physical_binding()
    interface = _binding(
        box_key=physical.box_key,
        source_id=physical.source_key[0],
        source_revision=physical.source_key[1],
        prepared_instance_id=physical.prepared_key[2],
    )

    linked = TheoremBackedSlowBoxBaseFieldBinding(
        interface=interface, physical_base=physical
    )

    assert linked.box_key == physical.box_key
    assert linked.source_key == physical.source_key
    assert linked.prepared_key == physical.prepared_key
    assert linked.physical_velocity_definition == PINNED_FINAL_SLOW_BASE_VELOCITY
    assert linked.theorem_physical_base_identity_linked is True
    assert linked.binary64_runtime_box_supported is False
    assert linked.status == "formal-structure"
    assert linked.actual_physical_base_values_materialized is False
    assert linked.base_jet_values_machine_certified is False
    assert linked.uniform_phase_estimates_proved is False
    assert linked.paper_exact_velocity_available is False


@pytest.mark.parametrize(
    ("field", "value", "match"),
    [
        ("box_key", (100000, (8, -3, 11)), "box key"),
        ("source_revision", "fixture-r2", "source identity"),
        ("prepared_instance_id", "different-prepared", "Prepared identity"),
    ],
)
def test_theorem_bridge_rejects_cross_wired_physical_identity(field, value, match):
    physical = _physical_binding()
    kwargs = dict(
        box_key=physical.box_key,
        source_id=physical.source_key[0],
        source_revision=physical.source_key[1],
        prepared_instance_id=physical.prepared_key[2],
    )
    kwargs[field] = value
    with pytest.raises(ValueError, match=match):
        TheoremBackedSlowBoxBaseFieldBinding(
            interface=_binding(**kwargs),
            physical_base=physical,
        )


def test_theorem_metadata_can_name_large_band_without_faking_binary64_runtime():
    binding = _binding(box_key=(100000, (7, -3, 11)))
    assert binding.box_key == (100000, (7, -3, 11))
    assert binding.binary64_runtime_box_supported is False

    with pytest.raises(ValueError, match="1<=ell<=1000"):
        SlowLabel(100000, (7, -3, 11), 1)


@pytest.mark.parametrize("producer", ["sampled", "fitted", "numeric-scan", "paper-exact"])
def test_non_interface_export_producers_are_rejected(producer):
    with pytest.raises(ValueError, match="formal-interface-export"):
        _binding(producer_kind=producer)


@pytest.mark.parametrize(
    "field",
    [
        "box_identity_certified",
        "sign_independence_certified",
        "normalized_chart_coordinates_certified",
        "provider_derivation_certified",
        "application_certified",
    ],
)
def test_missing_interface_certificate_fails_closed(field):
    with pytest.raises(ValueError, match="interface-export certified true"):
        _binding(**{field: False})


@pytest.mark.parametrize(
    "field",
    [
        "base_field_family_id",
        "normalized_base_field_id",
        "tangential_jet_provider_id",
        "evaluation_application_id",
        "source_id",
        "source_revision",
        "prepared_instance_id",
        "provenance",
    ],
)
def test_binding_identities_must_be_nonempty(field):
    with pytest.raises(ValueError, match=f"nonempty {field}"):
        _binding(**{field: " "})


@pytest.mark.parametrize(
    "box_key",
    [
        (0, (1, 2, 3)),
        (200, (1, 2)),
        (200, (1, True, 3)),
        ("200", (1, 2, 3)),
    ],
)
def test_malformed_box_keys_fail_closed(box_key):
    with pytest.raises(ValueError, match="box_key"):
        _binding(box_key)


def test_wrong_composition_types_are_rejected():
    rep = _active_representative()
    provider = Provider()
    with pytest.raises(TypeError, match="ActiveSlowRepresentative"):
        freeze_bound_label_phase_data(
            object(), provider, _binding(), rectangle_radius=0.01, u_star=0.25
        )
    with pytest.raises(TypeError, match="SlowBoxBaseFieldBindingWitness"):
        freeze_bound_label_phase_data(
            rep, provider, object(), rectangle_radius=0.01, u_star=0.25
        )

    good = freeze_bound_label_phase_data(
        rep, provider, _binding(), rectangle_radius=0.01, u_star=0.25
    )
    with pytest.raises(TypeError, match="FrozenLabelPhaseData"):
        BoundFrozenLabelPhaseData(object(), good.binding)
    with pytest.raises(TypeError, match="SlowBoxBaseFieldBindingWitness"):
        BoundFrozenLabelPhaseData(good.frozen, object())

    physical = _physical_binding()
    with pytest.raises(TypeError, match="SlowBoxBaseFieldBindingWitness"):
        TheoremBackedSlowBoxBaseFieldBinding(object(), physical)
    with pytest.raises(TypeError, match="LargeBandPhysicalBaseBinding"):
        TheoremBackedSlowBoxBaseFieldBinding(
            _binding(
                box_key=physical.box_key,
                source_id=physical.source_key[0],
                source_revision=physical.source_key[1],
                prepared_instance_id=physical.prepared_key[2],
            ),
            object(),
        )
