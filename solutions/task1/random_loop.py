"""
Implementation of random-based loop.

See solutions/run.py for usage.
"""
from pathlib import Path
from typing import List

from src.al_loop import Loop, LeadCompound
from src.compound_spaces import SmallZINC


class RandomLoop(Loop):
    """Samples random compounds from the ZINC database."""
    def __init__(self, base_dir: Path, user_token=None, target="DRD2"):
        self.space = SmallZINC()
        super().__init__(base_dir, user_token, target)

    def propose_candidates(self, n_candidates: int) -> List[LeadCompound]:
        smi = []
        while len(smi) != n_candidates:
            sampled_smi = self.space.try_sample()[0]
            if sampled_smi not in smi:
                smi.append(sampled_smi)
        return [
            LeadCompound(s, None, None) for s in smi
        ]