"""
Routes of the application
"""

import hashlib
import traceback

import numpy as np
import rdkit
from flask import request, render_template, jsonify
from flask_socketio import emit
from rdkit import Chem

from server.app import TOP_N, socketio, db, app
from server.models import Result, User
from src.eval import virtual_screen_TDC


def compute_md5(data):
    m = hashlib.md5()
    m.update(str(data).encode('utf-8'))
    return m.hexdigest()

@app.route("/results-checksum", methods=['GET'])
def results_checksum():
    results = Result.query.all()
    return jsonify({"checksum": compute_md5(results)})

@app.route('/score_compound', methods=['POST'])
def score_compounds():
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

        metrics = {}
        for k in [1, 10, 100]:
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