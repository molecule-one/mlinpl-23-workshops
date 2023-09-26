"""Standalone server to keep track of and serve the leaderboard"""
from flask import Flask
from flask_socketio import SocketIO
from flask_sqlalchemy import SQLAlchemy

TOP_N = 100
app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///results.db'
db = SQLAlchemy(app)
socketio = SocketIO(app)
# Define predefined call limits for each oracle_name. Adjust as needed.
call_limits = {
    "DRD2": 10000  # Example: 10 calls for DRD2 oracle. Add other oracles and limits as needed.
}

def sum_filter(value):
    return sum(value.values())

app.jinja_env.filters['sum'] = sum_filter