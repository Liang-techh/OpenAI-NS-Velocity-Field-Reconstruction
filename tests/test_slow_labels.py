import math

import numpy as np
import pytest

from openai_ns_reconstruction.charts import DyadicChart
from openai_ns_reconstruction.phase import TangentialBaseJet
from openai_ns_reconstruction.slow_labels import (
    ActiveSlowRepresentative,
    SlowBoxEnclosure,
    SlowLabel,
    freeze_label_phase_data,
)


def _active_representative(*, sigma=1, ratio=1.25):
    chart = DyadicChart(ell=200, h=0.005)
    R0 = 1.0
    Z0 = 0.2
    # Dimensionless form of q-z^2 q^(2h)=tau after q=Q*ratio.
    T0 = ratio - Z0**2 * ratio ** (2.0 * chart.h)
    label = SlowLabel(chart.ell, (7, -3, 11), sigma)
    mesh = chart.S_star ** -3
    enclosure = SlowBoxEnclosure(
        chart=chart,
        label=label,
        center=(R0, Z0, T0),
        half_width=(mesh / 2.0, mesh / 2.0, mesh / 2.0),
    )
    return ActiveSlowRepresentative(
        chart=chart,
        label=label,
        R0=R0,
        Z0=Z0,
        T0=T0,
        q_big=1.5 * chart.Q if ratio < 1.5 else 3.0 * chart.Q,
        X_a=0.35,
        X_b=0.45,
        enclosure=enclosure,
    )


def test_label_signs_duplicate_one_slow_box_and_mesh_scale_is_exact():
    rep = _active_representative(sigma=1)
    opposite = rep.label.opposite()
    assert opposite.sigma == -1
    assert opposite.box_key == rep.label.box_key
    assert rep.enclosure.mesh_size == pytest.approx(rep.chart.ell ** -6, rel=0.0, abs=0.0)

    with pytest.raises(ValueError, match="exceeds the paper"):
        SlowBoxEnclosure(
            chart=rep.chart,
            label=rep.label,
            center=(rep.R0, rep.Z0, rep.T0),
            half_width=(1.01 * rep.enclosure.mesh_size,) * 3,
        )


def test_active_representative_reconstructs_eq_6_8_dimensionless_shell():
    rep = _active_representative(ratio=1.25)
    assert rep.q_ratio == pytest.approx(1.25, rel=2e-13)
    assert rep.q / rep.chart.Q == pytest.approx(1.25, rel=2e-13)
    assert rep.X0 == pytest.approx(0.4, rel=2e-13)
    assert all(rep.active_hypotheses().values())
    assert np.array_equal(rep.slow_point, np.array([rep.R0, rep.Z0, rep.T0]))


def test_active_representative_fails_closed_outside_dyadic_band():
    chart = DyadicChart(ell=200, h=0.005)
    ratio = 2.1
    R0, Z0 = 1.0, 0.2
    T0 = ratio - Z0**2 * ratio ** (2.0 * chart.h)
    label = SlowLabel(chart.ell, (0, 0, 0), 1)
    mesh = chart.S_star ** -3
    enclosure = SlowBoxEnclosure(
        chart=chart,
        label=label,
        center=(R0, Z0, T0),
        half_width=(mesh / 2.0,) * 3,
    )
    with pytest.raises(ValueError, match="q_over_Q_band"):
        ActiveSlowRepresentative(
            chart=chart,
            label=label,
            R0=R0,
            Z0=Z0,
            T0=T0,
            q_big=3.0 * chart.Q,
            X_a=0.2,
            X_b=0.3,
            enclosure=enclosure,
        )


def test_freeze_label_phase_data_uses_only_the_fixed_base_representative():
    rep = _active_representative(sigma=-1)

    class Provider:
        def __init__(self):
            self.calls = []

        def tangential_jet(self, *, chart, R, Z, T):
            self.calls.append((chart.ell, R, Z, T))
            # At R0=1: g0=(R0*F_R,G_R)=(-3,0), which gives lambda0^2=2.
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

    provider = Provider()
    r0 = 0.01
    u_star = 0.25
    frozen = freeze_label_phase_data(
        rep, provider, rectangle_radius=r0, u_star=u_star
    )

    assert provider.calls == [(rep.chart.ell, rep.R0, rep.Z0, rep.T0)]
    assert frozen.frame.lambda0 == pytest.approx(math.sqrt(2.0))
    assert frozen.frame.c0 == pytest.approx(-math.sqrt(2.0) / 2.0)
    assert frozen.pulse_length == pytest.approx(2.0 * r0 / rep.chart.c_i)
    expected_k = math.ceil(rep.chart.epsilon ** -0.5)
    assert frozen.carrier_k == expected_k
    expected_B2 = math.sqrt(2.0) / (
        rep.chart.epsilon * expected_k**2 * (1.0 + u_star**2) ** 1.5
    )
    assert frozen.B_s**2 == pytest.approx(expected_B2)
    assert frozen.representative_phase.sigma == -1
    assert frozen.representative_phase.R0 == rep.R0
    assert np.allclose(frozen.representative_phase.g, np.array([-3.0, 0.0]))
    assert np.allclose(frozen.representative_phase.K, frozen.frame.K)


def test_freeze_label_phase_data_rejects_nonpaper_growth_jet_or_untyped_provider():
    rep = _active_representative()

    class BadGrowth:
        def tangential_jet(self, **kwargs):
            return TangentialBaseJet(F=1.0, G=0.0, F_R=1.0, G_R=0.0)

    with pytest.raises(ValueError, match="positive-growth"):
        freeze_label_phase_data(rep, BadGrowth(), rectangle_radius=0.01, u_star=0.2)

    class WrongType:
        def tangential_jet(self, **kwargs):
            return (1.0, 2.0)

    with pytest.raises(TypeError, match="TangentialBaseJet"):
        freeze_label_phase_data(rep, WrongType(), rectangle_radius=0.01, u_star=0.2)
