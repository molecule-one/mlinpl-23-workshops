"""
Implementation of the simplest virtual screening experiments.
"""
from typing import List, Dict, Literal, Tuple

import numpy as np
from tdc import Oracle
from src.utils import RdkitCanonicalSmiles

# TODO: Remove me
# import sklearn
# sys.modules["sklearn.svm.classes"] = sklearn.svm  # hotfix necessary due to TDC bug

def _get_TDC_oracle(oracle_name):
    return Oracle(name="DRD2")

def virtual_screen_TDC(
    compounds: List[RdkitCanonicalSmiles], oracle_name: str = "DRD2"
) -> List[float]:
    """
    Perform virtual screening in the space for compounds achieving high score according to a selected TDC oracle.
    """
    oracle = Oracle(name=oracle_name)
    return oracle(list(compounds))

def run_virtual_screening(compounds: List[RdkitCanonicalSmiles], experiment: Literal["DRD2"] ="DRD2") -> Tuple[Dict, List]:
    """Runs virtual screening for a list of spaces."""
    if experiment == "DRD2":
        fnc = virtual_screen_TDC
    else:
        raise NotImplementedError(f"Unknown experiment f{experiment}")

    scores = fnc(compounds)

    sorted_scores = sorted(scores)

    metrics = {}
    for k in [1, 10, 100]:
        metrics[f"top_{k}"] = np.mean(sorted_scores[-k:])

    return metrics, scores
