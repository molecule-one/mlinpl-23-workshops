"""Search at random in the space of compounds.

Usage: PORT=8000 python solutions/run.py -b 1000 -w random -t GSK3β
"""
import os
import shutil
from pathlib import Path
from typing import List

import argh
import numpy as np
from rich.console import Console
import matplotlib.pylab as plt

from src.   server_wrapper import FlaskAppClient
from src.al_loop import LeadCompound
from solutions.task1.random_loop import RandomLoop
from solutions.task2.mutate_loop import MutateLoop
from solutions.task4.ml_loop import MLLoop

console = Console()

BASEURL = os.environ.get("BASEURL", "http://127.0.0.1:8000")

def run(budget=1000, target="DRD2", purge=False, which="random", steps=10, user_token='test-0'):
    base_dir = Path("solution_{}_search_budget={}_target={}_user={}".format(which, budget, target, user_token))

    if purge:
        shutil.rmtree(base_dir, ignore_errors=True)

    if target == "GSK3β":
        client = None
    else:
        client = FlaskAppClient(BASEURL)

    if which == "random":
        loop = RandomLoop(base_dir=base_dir,
                          user_token=user_token,
                          target=target)
    elif which == "mutate":
        # NOTE: Exploration is key for mutate to work
        examples_per_step = budget // steps
        n_warmup_iterations = 300 // examples_per_step
        loop = MutateLoop(base_dir=base_dir,
                          n_warmup_iterations=n_warmup_iterations,
                          user_token=user_token,
                          target=target)
    elif which == "mutate_low_exp":
        # NOTE: Exploration is key for mutate to work
        loop = MutateLoop(base_dir=base_dir,
                          n_warmup_iterations=1,
                          user_token=user_token,
                          target=target)
    elif which == "ml":
        examples_per_step = budget // steps
        n_warmup_iterations = 300 // examples_per_step
        base_loop = MutateLoop(base_dir=base_dir, n_warmup_iterations=n_warmup_iterations, user_token=user_token, target=target)
        loop = MLLoop(base_dir=base_dir, n_warmup_iterations=n_warmup_iterations*2, base_loop=base_loop, user_token=user_token, target=target)
    else:
        raise ValueError(f"Unknown which={which}")

    if loop.n_iterations > 0:
        raise ValueError("Already run. Please remove the folder solution_random_search to run again.")

    metrics = []
    all_result: List[LeadCompound] = []
    budget_per_step = budget // steps
    assert budget % steps == 0 # for simplicity
    for step in range(steps):
        console.print(f"[red]Step {step}[/red]")
        candidates = loop.propose_candidates(budget_per_step)
        loop.test_in_lab_and_save(candidates, client=client)
        result: List[LeadCompound] = loop.load(iteration_id=step)
        all_result += result
        loop.generate_visualization(iteration_id=step)
        all_result_sorted = sorted(all_result, key=lambda x: x.activity, reverse=True)
        metrics.append({"top10": np.mean([x.activity for x in all_result_sorted[:10]]),
                        "top10_synth": np.mean([x.synth_score for x in all_result_sorted[:10]])})

    # plot metrics using matplotlib and save
    for k in metrics[0]:
        plt.plot([i*budget_per_step for i in range(len(metrics))],
                 [m[k] for m in metrics], linewidth=4)
        plt.grid()
        plt.yticks(np.arange(0, 1.1, 0.1))
        plt.xlabel('N_compounds')
        plt.ylabel(k)
        plt.savefig(base_dir / f"{k}.png")
        plt.close()


if __name__ == "__main__":
    argh.dispatch_command(run)