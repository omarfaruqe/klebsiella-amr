import math

import pytest

from kamr.analysis.predictive_values import predictive_values


def test_perfect_test_is_perfectly_predictive():
    ppv, npv = predictive_values(1.0, 1.0, 0.1)
    assert ppv == 1.0
    assert npv == 1.0


def test_ppv_falls_as_prevalence_falls():
    # the point of the limitation: discrimination is fixed, predictive value is not
    high = predictive_values(0.567, 0.952, 0.24)[0]
    low = predictive_values(0.567, 0.952, 0.02)[0]
    assert high > low
    assert low < 0.25


def test_npv_rises_as_prevalence_falls():
    assert predictive_values(0.567, 0.952, 0.02)[1] > predictive_values(0.567, 0.952, 0.24)[1]


def test_matches_bayes_by_hand():
    ppv, npv = predictive_values(0.6, 0.9, 0.1)
    assert math.isclose(ppv, 0.06 / (0.06 + 0.09))
    assert math.isclose(npv, 0.81 / (0.81 + 0.04))


def test_uninformative_test_returns_prevalence_as_ppv():
    # a rule that calls everything resistant has PPV equal to prevalence
    ppv, _ = predictive_values(1.0, 0.0, 0.3)
    assert math.isclose(ppv, 0.3)


@pytest.mark.parametrize("bad", [-0.1, 1.1])
def test_rejects_out_of_range_inputs(bad):
    with pytest.raises(ValueError):
        predictive_values(bad, 0.9, 0.1)
    with pytest.raises(ValueError):
        predictive_values(0.9, 0.9, bad)
