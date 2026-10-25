import numpy as np
import pytest

from kamr.benchmark.compare import bh_qvalues, bootstrap_pvalue


def test_pvalue_is_one_when_the_distribution_straddles_zero_evenly():
    diffs = np.array([-2.0, -1.0, 1.0, 2.0])
    assert bootstrap_pvalue(diffs) == 1.0


def test_pvalue_hits_the_resolution_floor_when_no_replicate_crosses_zero():
    diffs = np.ones(1000)
    # with 1000 resamples nothing below 1/1000 is distinguishable, and it must not be reported
    # as zero, which would imply certainty the bootstrap cannot supply
    assert bootstrap_pvalue(diffs) == pytest.approx(1 / 1000)


def test_pvalue_is_symmetric_in_sign():
    diffs = np.concatenate([np.ones(950), -np.ones(50)])
    assert bootstrap_pvalue(diffs) == pytest.approx(bootstrap_pvalue(-diffs))


def test_pvalue_of_empty_distribution_is_one():
    assert bootstrap_pvalue(np.array([])) == 1.0


def test_bh_leaves_a_single_test_unchanged():
    assert bh_qvalues(np.array([0.03])) == pytest.approx([0.03])


def test_bh_matches_worked_example():
    # classic BH: p * m / rank, then enforced monotone from the top
    p = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
    expected = [0.05, 0.05, 0.05, 0.05, 0.05]
    assert bh_qvalues(p) == pytest.approx(expected)


def test_bh_preserves_input_order():
    p = np.array([0.5, 0.001, 0.2])
    q = bh_qvalues(p)
    assert np.argmin(q) == 1


def test_bh_is_monotone_and_never_exceeds_one():
    rng = np.random.default_rng(0)
    p = rng.random(200)
    q = bh_qvalues(p)
    assert q.max() <= 1.0
    # sorted by p, q must be non-decreasing
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)


def test_bh_is_less_conservative_than_bonferroni():
    p = np.array([0.001, 0.01, 0.02, 0.04])
    q = bh_qvalues(p)
    assert np.all(q <= p * len(p) + 1e-12)
