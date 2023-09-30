"""
Routes of the application
"""

import hashlib
import traceback
from typing import List
from uuid import uuid4

import numpy as np
import rdkit
from flask import request, render_template, jsonify
from flask_socketio import emit
from rdkit import Chem

from sas_score import compute_ertl_score
from server.app import TOP_N, socketio, db, app, MASTER_KEY, call_limits
from server.models import Result, User, Token
from src.eval import virtual_screen_TDC

from rich.console import Console
console = Console()

def _compute_md5(data):
    m = hashlib.md5()
    m.update(str(data).encode('utf-8'))
    return m.hexdigest()

def _validate_smiles(cls, candidates: List[str]):
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

def _evaluate_synthesizability(cls, candidates: List[str]) -> List[float]:
    cls._validate_smiles([c for c in candidates])
    return [compute_ertl_score(c) for c in candidates]

@app.route("/results-checksum", methods=['GET'])
def results_checksum():
    results = Result.query.all()
    return jsonify({"checksum": _compute_md5(results)})

@app.route('/score_compound', methods=['POST'])
def score_compound():
    try:
        compound = request.json.get('compound')

        try:
            mol = rdkit.Chem.MolFromSmiles(compound)
            if mol is None:
                return jsonify({"error": f"Failed to parse SMILES {compound} (rdkit.Chem.MolFromSmiles(smi) returns None)."}), 500
        except:
            # Get the traceback details and return it along with the error message
            tb = traceback.format_exc()
            return jsonify({"error": f"Failed to parse SMILES {compound} (rdkit.Chem.MolFromSmiles(smi) throws an error).", "traceback": tb}), 500

        oracle_name = request.json.get('oracle_name', "DRD2")
        token = request.json.get('token')

        user = User.query.get(token)
        if not user:
            user = User(id=token, oracle_calls={})
            db.session.add(user)
        if oracle_name not in user.oracle_calls:
            user.oracle_calls[oracle_name] = 0

        if user.oracle_calls[oracle_name] >= call_limits.get(oracle_name, float('inf')):
            return jsonify({"error": f"Call limit reached for oracle: {oracle_name}"}), 403

        score = virtual_screen_TDC([compound], oracle_name)[0]
        user.oracle_calls[oracle_name] += 1
        db.session.commit()

        return jsonify({"score": score}), 200

    except Exception as e:
        # Get the traceback details and return it along with the error message
        tb = traceback.format_exc()
        return jsonify({"error": str(e), "traceback": tb}), 500

@app.route('/leaderboard')
def index():
    results = Result.query.all()
    sorted_results = sorted(results, key=lambda x: -sum(x.metrics.values()))
    return render_template('index.html', results=sorted_results)

@app.route("/evaluate_and_add_result", methods=['POST'])
def evaluate_and_add_result():
    try:
        token = request.json.get('token')

        # Check if the token is valid
        if not Token.check_valid_token(token):
            return jsonify({"status": "failure", "message": "Invalid token"}), 403

        compounds = request.json.get('compounds').split(",")
        if len(compounds) > TOP_N:
            compounds = np.random.RandomState(777).choice(compounds, TOP_N)
        oracle_name = request.json.get('oracle_name')

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
        scores = virtual_screen_TDC(compounds, oracle_name)
        scores = sorted(scores)

        # note: duplication of code of add_result
        metrics = {}
        for k in [10]:
            metrics[f"{oracle_name}_top_{k}"] = np.mean(scores[-k:])

        result = Result.query.get(token)

        if not result:
            result = Result(id=token, metrics=metrics)
            db.session.add(result)
        else:
            result.metrics.update(metrics)

        db.session.commit()
        return jsonify({"status": "success", "metrics": metrics}), 200

    except Exception as e:
        # Get the traceback details and return it along with the error message
        tb = traceback.format_exc()
        return jsonify({"error": str(e), "traceback": tb}), 500

@app.route('/add_result', methods=['POST'])
def add_result():
    token = request.json.get('token')
    metrics = request.json.get('metrics')

    # Check if the token is valid
    if not Token.check_valid_token(token):
        return jsonify({"status": "failure", "message": "Invalid token"}), 403

    result = Result.query.get(token)

    if not result:
        result = Result(id=token, metrics=metrics)
        db.session.add(result)
    else:
        result.metrics = metrics

    db.session.commit()

    return jsonify({"status": "success"}), 200

@app.route('/reset', methods=['POST'])
def reset():
    db.session.query(Result).delete()
    db.session.commit()

    socketio.emit('update', {})
    return jsonify({"status": "reset completed"}), 200

@socketio.on('connect')
def test_connect():
    emit('update', {})


@app.route('/generate_tokens', methods=['POST'])
def generate_tokens():
    # This is an example for demonstration. In a real-world scenario, you might want to secure this endpoint.
    master_key = request.json.get('master_key')
    if master_key != MASTER_KEY:  # Replace with your secret key
        return jsonify({"status": "failure", "message": "Invalid master key"}), 403

    tokens = []
    for _ in range(50):
        new_token = str(uuid4())  # generate a unique token
        tokens.append(new_token)
        db.session.add(Token(token=new_token))

    db.session.commit()

    return jsonify({"status": "success", "tokens": tokens}), 200


@app.route('/list_tokens', methods=['GET'])
def list_tokens():
    all_tokens = Token.query.all()
    token_list = [t.token for t in all_tokens]
    return jsonify(token_list), 200
