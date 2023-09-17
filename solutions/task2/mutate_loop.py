"""
Implementation of AL algorithm that mutates top compounds from the previous iterations.

See solutions/run.py for usage.
"""
from pathlib import Path
from typing import List

import numpy as np

from mlinpl_workshops.src.al_loop import Loop, LeadCompound
from mlinpl_workshops.src.compound_spaces import SmallZINC

from selfies import encoder, decoder

from rich.console import Console

from src.mutate import mutate_selfie

console = Console()

class MutateLoop(Loop):
    """Implementation of AL algorithm that mutates top compounds from the previous iterations."""
    def __init__(self, base_dir: Path, n_warmup_iterations: int = 1, mutate_top_k: int = 10, user_token=None, target="DRD2"):
        self.space = SmallZINC()
        self.n_warmup_iterations = n_warmup_iterations
        self.mutate_top_k = mutate_top_k
        super().__init__(base_dir, user_token, target)

    def _propose_random(self, n_candidates: int) -> List[LeadCompound]:
        smi = [self.space.try_sample()[0] for _ in range(n_candidates)]
        return [
            LeadCompound(s, None, None) for s in smi
        ]

    def propose_candidates(self, n_candidates: int) -> List[LeadCompound]:
        previous_results: List[LeadCompound] = self.load()

        if n_candidates < self.mutate_top_k:
            raise ValueError(f"n_candidates must be at least mutate_top_k ({self.mutate_top_k}).")

        if n_candidates == 0:
            return []

        if self.n_iterations < self.n_warmup_iterations:
            return self._propose_random(n_candidates)

        topK_ids = np.argsort([c.activity for c in previous_results])[-self.mutate_top_k:]
        console.log("Mutating top compounds:")
        for i in topK_ids:
            console.log(previous_results[i])
        selfies = [encoder(previous_results[i].smiles) for i in topK_ids]
        m = n_candidates // self.mutate_top_k + 1
        new_compounds = []
        for i in range(self.mutate_top_k):
            new_smiles = []
            m_target = min(m, n_candidates - len(new_compounds)) # last batch might be smaller
            while len(new_smiles) != m_target:
                if len(new_compounds) == n_candidates:
                    break
                new_smi = decoder(mutate_selfie(selfies[i], max_molecules_len=100)[0])
                if new_smi not in new_smiles and new_smi not in new_compounds:
                    new_smiles.append(new_smi)
            new_compounds += new_smiles
        assert len(new_compounds) == n_candidates
        return [LeadCompound(smiles=c) for c in new_compounds]