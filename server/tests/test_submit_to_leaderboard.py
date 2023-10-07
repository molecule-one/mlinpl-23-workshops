"""
Tests for the server.

Before running them, the server database should be restarted.

Run as: py.test server/tests -s
"""

import shutil
from pathlib import Path
from typing import List

import numpy as np

from src.al_loop import LeadCompound
from server.app import WORKSHOP_ORACLES
from solutions.task1.random_loop import RandomLoop
from src.server_wrapper import FlaskAppClient

from rich.console import Console

console = Console()

TEST_TOKEN_PREFIX = 'test-' # test-0, test-1, ...

def test_submitting_compounds_to_workshop_oracles():
    """Submits three simple molecules to the server using token test-0, to all workshop oracles."""
    client = FlaskAppClient()
    token = TEST_TOKEN_PREFIX + '0'

    # Example for scoring compounds

    for oracle in WORKSHOP_ORACLES:
        compounds = ["CCCCCCCCC", "CCCCCCCC", "CCCCCC=O"]
        response = client.score_compounds_and_update_leaderboard(compounds, oracle, token)
        print(response)
        assert "metrics" in response
        assert "compound_scores" in response
        assert "compound_sas_scores" in response


def test_random_exploration_gets_reasonable_score():
    """Simple random exploration of ZINC. Should get above >0.5 score on each oracle."""
    base_dir = Path("tmp")
    for protein in ['DRD2_server', 'JNK3', 'GSK3β']:
        shutil.rmtree(base_dir, ignore_errors=True)
        loop = RandomLoop(base_dir=base_dir,
                              user_token='test-1',
                              target=protein)
        all_result: List[LeadCompound] = []
        budget_per_step = 100
        for step in range(10):
            console.print(f"[red]Step {step}[/red]")
            candidates = loop.propose_candidates(budget_per_step)
            loop.test_in_lab_and_save(candidates)
            result: List[LeadCompound] = loop.load(iteration_id=step)
            all_result += result
            all_result_sorted = sorted(all_result, key=lambda x: x.activity, reverse=True)
            metrics = {"top10": np.mean([x.activity for x in all_result_sorted[:10]]),
                            "top10_synth": np.mean([x.synth_score for x in all_result_sorted[:10]])}
            console.log(metrics)

        assert metrics['top10'] > 0.1, "Random search should identify reasonable compounds"

test_random_exploration_gets_reasonable_score()