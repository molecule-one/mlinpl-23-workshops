"""
Solution to task 3

It loads previous results using Loop base class and simply recomputes Top10.

Usage: python solutions/task3/compare_random_mutate.py
"""

from pathlib import Path
from typing import List

import numpy as np
from matplotlib import pyplot as plt

from solutions.task1.random_loop import RandomLoop
from solutions.task2.mutate_loop import MutateLoop
from rich.console import Console

from src.al_loop import LeadCompound

console = Console()

if __name__ == "__main__":
    random_base_dir = Path("solution_{}_search_budget={}".format('random', 1000))
    random_loop = RandomLoop(base_dir=random_base_dir)
    mutate_base_dir = Path("solution_{}_search_budget={}".format('mutate', 1000))
    mutate_loop = MutateLoop(base_dir=mutate_base_dir)
    budget_per_step = 100

    for name, loop in [("random", random_loop), ("mutate", mutate_loop)]:
        all_result: List[LeadCompound] = []
        metrics = []
        for step in range(10):
            result: List[LeadCompound] = loop.load(iteration_id=step)
            all_result += result
            loop.generate_visualization(iteration_id=step)
            all_result_sorted = sorted(all_result, key=lambda x: x.activity, reverse=True)
            metrics.append({"top10": np.mean([x.activity for x in all_result_sorted[:10]]),
                            "top10_synth": np.mean([x.synth_score for x in all_result_sorted[:10]])})

        plt.grid(True)
        plt.plot([i * budget_per_step for i in range(len(metrics))],
                 [m["top10"] for m in metrics], linewidth=4, label=name)
        plt.xlabel('N_compounds')
        plt.ylabel('Top 10')
    plt.savefig(f"mutate_vs_random.png")
    plt.close()