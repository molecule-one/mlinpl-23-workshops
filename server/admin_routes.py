"""
Routes of the application that are only visible to admin
"""

import hashlib
from uuid import uuid4

from flask import request, jsonify
from flask_socketio import emit
from rich.console import Console

from server.app import socketio, db, app, MASTER_KEY
from server.models import Result, Token, User

console = Console()

def _compute_md5(data):
    m = hashlib.md5()
    m.update(str(data).encode('utf-8'))
    return m.hexdigest()


@app.route("/results-checksum", methods=['GET'])
def results_checksum():
    results = Result.query.all()
    return jsonify({"checksum": _compute_md5(results)})

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
    db.session.query(User).delete()
    db.session.query(Result).delete()
    db.session.query(Token).delete()
    db.session.commit()

    # initiate database with 50 tokens of kind test-0, ... test-10
    tokens = []
    for i in range(50):
        new_token = "test-" + str(i)
        tokens.append(new_token)
        db.session.add(Token(token=new_token))

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
