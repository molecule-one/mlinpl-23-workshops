"""
Solution to task 3

It loads previous results using Loop base class and simply recomputes Top10.

Usage:
    python solutions/task3/compare.py -e solution_random_search_budget=1000,solution_mutate_search_budget=1000 -n random,mutate
"""
from functools import partial
from pathlib import Path
from typing import List

import argh
import numpy as np
from matplotlib import pyplot as plt
from more_itertools import zip_equal

from solutions.task1.random_loop import RandomLoop
from solutions.task2.mutate_loop import MutateLoop
from rich.console import Console

from src.al_loop import LeadCompound

from solutions.task4.ml_loop import MLLoop

console = Console()

def _get_cls(exp_name):
    # hacky function to determine cls of the Loop from exp name
    if "random" in exp_name:
        return RandomLoop
    elif "mutate" in exp_name:
        return MutateLoop
    elif "ml" in exp_name:
        return partial(MLLoop, base_loop=None) # hack just to initialize the class
    else:
        raise ValueError(f"Unknown exp_name={exp_name}")

def compare(exps="solution_{}_search_budget={}".format('random', 1000) + "," +  "solution_{}_search_budget={}".format('mutate', 1000),
            names="random,mutate"):
    exps = exps.split(",")
    names = names.split(",")
    assert len(exps) == len(names)
    loops = [_get_cls(exp)(base_dir=Path(exp)) for exp in exps]
    budget_per_step = 100 # TODO: might not always be the case

    for name, loop in zip_equal(names, loops):
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
    plt.title("Top10 comparison")
    plt.legend()
    plt.savefig("_".join(exps) + ".png")
    plt.close()

if __name__ == '__main__':
    argh.dispatch_command(compare)