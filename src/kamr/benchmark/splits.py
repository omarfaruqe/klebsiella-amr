"""Train/test split schemes for the benchmark: random, study-aware, lineage-aware, temporal.

Each scheme yields (train_index, test_index) pairs of positional indices into a label table.
Grouped schemes give every genome with a missing group value its own group, so missing metadata
never merges unrelated genomes. The temporal scheme uses dated genomes only.
"""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold, StratifiedKFold

SCHEMES = ("random", "study", "lineage", "temporal")
GROUP_COLUMN = {"study": "bioproject_accession", "lineage": "st"}


def groups_for(lab: pd.DataFrame, column: str) -> np.ndarray:
    values = lab[column].astype("object")
    fallback = "missing:" + lab["genome_id"].astype(str)
    return (
        values.where(values.notna() & (values.astype(str) != ""), fallback).astype(str).to_numpy()
    )


def temporal_split(
    lab: pd.DataFrame, test_share: float = 0.2
) -> tuple[np.ndarray, np.ndarray, int]:
    """Train on genomes up to the cut-off year, test on later ones (about test_share of dated)."""
    years = pd.to_numeric(lab["collection_year"], errors="coerce")
    dated = np.flatnonzero(years.notna().to_numpy())
    if len(dated) == 0:
        return np.array([], int), np.array([], int), -1
    y = years.iloc[dated].astype(int)
    cutoff = int(np.quantile(y, 1 - test_share))
    if (y > cutoff).sum() == 0:
        cutoff -= 1
    train = dated[(y <= cutoff).to_numpy()]
    test = dated[(y > cutoff).to_numpy()]
    return train, test, cutoff


def folds(lab: pd.DataFrame, scheme: str, n_splits: int = 5, seed: int = 20261024) -> Iterator:
    y = lab["y"].to_numpy()
    if scheme == "random":
        yield from StratifiedKFold(n_splits, shuffle=True, random_state=seed).split(lab, y)
    elif scheme in GROUP_COLUMN:
        g = groups_for(lab, GROUP_COLUMN[scheme])
        k = min(n_splits, len(set(g)))
        yield from GroupKFold(k).split(lab, y, g)
    elif scheme == "temporal":
        train, test, _ = temporal_split(lab)
        if len(train) and len(test):
            yield train, test
    else:
        raise ValueError(f"unknown scheme {scheme}")
