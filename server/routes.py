"""
Routes of the application:

* /leaderboard: shows up the leaderboard
* /score_compounds_and_update_leaderboard: scores provided compounds and updates the leaderboard
"""

import traceback
from typing import List

import numpy as np
import rdkit
from flask import request, render_template, jsonify
from more_itertools import zip_equal
from rdkit import Chem
from rich.console import Console

from src.sas_score import compute_ertl_score
from server.app import TOP_N, db, app, call_limits, SAS_THRESHOLD, WORKSHOP_ORACLES
from server.models import Result, User, Token
from src.eval import virtual_screen_TDC

console = Console()


def _validate_smiles(candidates: List[str]):
    """Helper function to check if the SMILES are valid"""
    for s in candidates:
        if not isinstance(s, str):
            raise ValueError("SMILES must be a string.")
        if len(s) == 0:
            raise ValueError("SMILES cannot be empty.")

        try:
            mol = rdkit.Chem.MolFromSmiles(s)
            if mol is None:
                raise ValueError("Invalid SMILES")
        except Exception as e:
            console.print_exception(show_locals=True)
            raise ValueError(f"Failed to parse SMILES using rdkit: {e}")

def _evaluate_synthesizability(candidates: List[str]) -> List[float]:
    _validate_smiles([c for c in candidates])
    return [compute_ertl_score(c) for c in candidates]

@app.route('/leaderboard')
def index():
    results = Result.query.all()
    sorted_results = sorted(results, key=lambda x: -sum(x.metrics.values()))
    return render_template('index.html', results=sorted_results)

@app.route("/score_compounds_and_update_leaderboard", methods=['POST'])
def score_compounds_and_update_leaderboard():
    """
    Scores compounds and updates leaderboard with the running top N score

    Compounds that are too hard to synthesize have returned score -1
    """
    try:
        token = request.json.get('token')
        oracle_name = request.json.get('oracle_name', "_DRD2")
        oracle_name = oracle_name.replace("_server", "")

        if oracle_name not in WORKSHOP_ORACLES:
            return jsonify({"status": "failure", "message": f"Expected oracle in {WORKSHOP_ORACLES}"}), 403

        # Check if the token is valid
        if not Token.check_valid_token(token):
            return jsonify({"status": "failure", "message": "Invalid token"}), 403

        user = User.query.get(token)
        if not user:
            user = User(id=token, oracle_calls={}, compound_scores={}, compound_sas_scores={})
            db.session.add(user)
        if oracle_name not in user.oracle_calls:
            user.oracle_calls[oracle_name] = 0

        n_remaining_calls = call_limits.get(oracle_name, float('inf')) - user.oracle_calls[oracle_name]

        if n_remaining_calls <= 0:
            return jsonify({"error": f"Call limit reached for oracle: {oracle_name}"}), 403

        compounds = request.json.get('compounds')
        if compounds is None:
            return jsonify({"error": "Missing 'compounds' field in the request."}), 500
        compounds = compounds.split(",")
        if len(compounds) > n_remaining_calls:
            compounds = np.random.RandomState(777).choice(compounds, n_remaining_calls)

        # update the limit
        user.oracle_calls[oracle_name] += len(compounds)

        for compound in compounds:
            try:
                mol = rdkit.Chem.MolFromSmiles(compound)
                if mol is None:
                    return jsonify(
                        {"error": f"Failed to parse SMILES {compound} (rdkit.Chem.MolFromSmiles(smi) returns None)."}), 500
            except:
                # Get the traceback details and return it along with the error message
                tb = traceback.format_exc()
                return jsonify({"error": f"Failed to parse SMILES {compound} (rdkit.Chem.MolFromSmiles(smi) throws an error).",
                                "traceback": tb}), 500
        # HACK: replaces "_server" which is special sequence to differntiate DRD2 from DRD2_server
        sas_scores = _evaluate_synthesizability(compounds)
        vs_scores = virtual_screen_TDC(compounds, oracle_name)
        scores = [vs_score if sas_score <= SAS_THRESHOLD else -1 for vs_score, sas_score in zip_equal(vs_scores, sas_scores)]

        if user.compound_scores is None:
            user.compound_scores = {}

        if oracle_name not in user.compound_scores:
            user.compound_scores[oracle_name] = []
            user.compound_sas_scores[oracle_name] = []

        user.compound_scores[oracle_name] += scores
        user.compound_sas_scores[oracle_name] += sas_scores

        # note: duplication of code of add_result
        metrics = {}
        for k in [TOP_N]:
            for oracle_name in WORKSHOP_ORACLES:
                if oracle_name in user.compound_scores:
                    top_ids = np.argsort(user.compound_scores[oracle_name])[-k:]
                    metrics[f"{oracle_name}_top_{k}"] = np.mean([user.compound_scores[oracle_name][i] for i in top_ids])
                else:
                    metrics[f"{oracle_name}_top_{k}"] = 0.0

        result = Result.query.get(token)

        if not result:
            result = Result(id=token, metrics=metrics)
            db.session.add(result)
        else:
            result.metrics.update(metrics)

        db.session.commit()
        return jsonify({"status": "success", "metrics": metrics, "compound_scores": scores, "compound_sas_scores": sas_scores}), 200

    except Exception as e:
        # Get the traceback details and return it along with the error message
        tb = traceback.format_exc()
        console.log(tb)
        return jsonify({"error": str(e), "traceback": tb}), 500
