import math

import pytest

from openai_ns_reconstruction.charts import DyadicChart
from openai_ns_reconstruction.slow_labels import SlowLabel
from openai_ns_reconstruction.slow_partition import (
    ProductSlowCutoff,
    active_dyadic_indices,
    active_product_indices,
    dyadic_cutoff,
    normalized_translate,
    product_partition_square_sum,
    slow_label_cutoff,
)


def test_normalized_integer_translates_form_local_squared_partition():
    for x in (-3.2, -1.5, -0.1, 0.0, 0.37, 0.5, 0.999, 7.25):
        i0 = math.floor(x)
        weights = [normalized_translate(x, i0), normalized_translate(x, i0 + 1)]
        assert sum(w * w for w in weights) == pytest.approx(1.0, rel=0.0, abs=8e-16)
        assert normalized_translate(x, i0 - 1) == 0.0
        assert normalized_translate(x, i0 + 2) == 0.0


def test_dyadic_partition_has_exact_paper_support_ratio_and_square_sum():
    ell = 40
    Q = math.ldexp(1.0, -ell)
    for ratio in (0.51, 0.75, 1.0, 1.3, 1.99):
        q = Q * ratio
        active = active_dyadic_indices(q)
        total = sum(dyadic_cutoff(j, q) ** 2 for j in active)
        assert total == pytest.approx(1.0, rel=0.0, abs=8e-16)

    assert dyadic_cutoff(ell, 0.5 * Q) == 0.0
    assert dyadic_cutoff(ell, 2.0 * Q) == 0.0
    assert dyadic_cutoff(ell, 0.49 * Q) == 0.0
    assert dyadic_cutoff(ell, 2.01 * Q) == 0.0


def test_product_partition_is_finite_local_and_sums_to_one():
    chart = DyadicChart(ell=10, h=0.005)
    hmesh = chart.S_star ** -3
    point = (
        (1_000_000 + 0.23) * hmesh,
        (-250_000 + 0.61) * hmesh,
        (500_000 + 0.44) * hmesh,
    )
    active = active_product_indices(chart, *point)
    assert len(active) == 8
    assert product_partition_square_sum(chart, *point) == pytest.approx(
        1.0, rel=0.0, abs=2e-15
    )
    assert all(ProductSlowCutoff(chart, a).support_contains(*point) for a in active)


def test_actual_product_cutoff_bridges_to_slow_box_enclosure_and_duplicates_sign():
    chart = DyadicChart(ell=10, h=0.005)
    a = (1_000_000, -250_000, 500_000)
    product = ProductSlowCutoff(chart, a)
    plus = SlowLabel(chart.ell, a, 1)
    minus = plus.opposite()
    enclosure = product.enclosure(plus)

    assert enclosure.center == product.center
    assert enclosure.half_width == (product.mesh_size,) * 3
    assert enclosure.contains(*product.center)

    Q = chart.Q
    q = 1.2 * Q
    value_plus = slow_label_cutoff(chart, plus, q, *product.center)
    value_minus = slow_label_cutoff(chart, minus, q, *product.center)
    assert value_plus == value_minus
    assert value_plus > 0.0


def test_partition_interfaces_fail_closed_on_invalid_inputs_or_mismatched_label():
    chart = DyadicChart(ell=10, h=0.005)
    with pytest.raises(ValueError, match="q must be positive"):
        dyadic_cutoff(chart.ell, 0.0)
    with pytest.raises(ValueError, match="integer"):
        normalized_translate(0.2, True)

    product = ProductSlowCutoff(chart, (1_000_000, 0, 0))
    wrong = SlowLabel(chart.ell + 1, product.a, 1)
    with pytest.raises(ValueError, match="label must identify"):
        product.enclosure(wrong)
    with pytest.raises(ValueError, match="same dyadic band"):
        slow_label_cutoff(chart, wrong, chart.Q, *product.center)
