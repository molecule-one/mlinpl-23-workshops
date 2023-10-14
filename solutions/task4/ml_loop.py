"""
Implementation of AL algorithm that uses ML algorithms to guide the search.
"""
from pathlib import Path
from typing import List

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

from src.al_loop import Loop, LeadCompound
from rich.console import Console

console = Console()

class MLLoop(Loop):
    """
    Implementation of AL algorithm that uses ML algorithms to guide the search.

    It trains RF internally.

    The algorithm:
        1. Retrain on all compounds so far
        2. Propose 10 x n_candidates compounds using base loop
        3. Evaluate them using ML models.
        4. Return n_candidates with the highest predicted activity.
    """
    def __init__(self, base_dir: Path, base_loop: Loop, n_warmup_iterations: int=1, user_token=None, target="DRD2"):
        self.base_loop = base_loop
        self.n_warmup_iterations = n_warmup_iterations
        super().__init__(base_dir, user_token, target)

    def _featurize(self, smi: List[str]) -> np.ndarray:
        mols = [Chem.MolFromSmiles(s) for s in smi]
        fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in mols]
        X = []
        for fp in fps:
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(fp, arr)
            X.append(arr)
        return np.array(X)

    def _train_model(self, previous_results: List[LeadCompound]):
        """Trains models and assigns to ._model variable."""
        previous_results = [c for c in previous_results if c.activity != -1]
        if len(previous_results) == 0:
            raise ValueError("No previous results to train on (excluded activity = -1). Perhaps your "
                             "base loop proposes nonsynthetizable compounds?")
        console.log(f"Retraining model on {len(previous_results)} compounds.")
        smi = [c.smiles for c in previous_results]
        X = self._featurize(smi)
        y = np.array([c.activity for c in previous_results])
        y = y > np.median(y) # convert to binary. loses information naturally.

        # split using sklearn
        X_temp, X_test, y_temp, y_test, smi_temp, smi_test = \
            train_test_split(X, y, smi, test_size=0.3, random_state=42)
        X_train, X_valid, y_train, y_valid, smi_train, smi_valid = \
            train_test_split(X_temp, y_temp, smi_temp, test_size=0.3, random_state=42)

        console.log(f"Training set size: {len(X_train)}")
        console.log(f"Validation set size: {len(X_valid)}")
        console.log(f"Test set size: {len(X_test)}")
        console.log(f"Training set activity mean: {np.mean(y_train)}")
        console.log(f"Proceeding to training Random Forest")
        rf = RandomForestClassifier(n_estimators=100, random_state=42)
        rf.fit(X_train, y_train)

        self._model = rf

    def _select_top_N(self, candidates: List[LeadCompound], n_select: int) -> List[LeadCompound]:
        """Ranks candidates by their predicted activity."""
        X_test = self._featurize([c.smiles for c in candidates])
        y_pred = self._model.predict_proba(X_test)
        if y_pred.ndim == 2:
            y_pred = y_pred[:, -1]
        return [c for _, c in sorted(zip(y_pred, candidates), reverse=True, key=lambda a: a[0])][:n_select]

    def propose_candidates(self, n_candidates: int) -> List[LeadCompound]:
        if self.n_iterations < self.n_warmup_iterations:
            return self.base_loop.propose_candidates(n_candidates)

        previous_results: List[LeadCompound] = self.load()
        self._train_model(previous_results)
        return self._select_top_N(self.base_loop.propose_candidates(10 * n_candidates), n_candidates)