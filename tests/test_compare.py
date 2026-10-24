import numpy as np
import pandas as pd

from kamr.benchmark.compare import balanced_accuracy, paired_bootstrap


def test_balanced_accuracy():
    y = np.array([1, 1, 0, 0])
    assert balanced_accuracy(y, np.array([1, 1, 0, 0])) == 1.0
    assert balanced_accuracy(y, np.array([1, 0, 0, 1])) == 0.5
    assert np.isnan(balanced_accuracy(np.ones(4, int), np.ones(4, int)))


def test_paired_bootstrap_detects_real_and_absent_differences():
    rng = np.random.default_rng(0)
    n = 400
    y = np.repeat([1, 0], n // 2)
    strong = y.copy()  # perfect predictor
    weak = y.copy()
    weak[:80] = 0  # misses 40% of positives
    wide = pd.DataFrame({"y": y, "a": strong, "b": weak, "c": strong})
    better = paired_bootstrap(wide, "a", "b", 400, rng)
    assert better["diff"] > 0.15 and better["excludes_zero"]
    same = paired_bootstrap(wide, "a", "c", 400, rng)
    assert same["diff"] == 0 and not same["excludes_zero"]
