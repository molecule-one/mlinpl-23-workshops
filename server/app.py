"""Standalone server to keep track of and serve the leaderboard"""
from flask import Flask
from flask_socketio import SocketIO
from flask_sqlalchemy import SQLAlchemy

# CONFIGURATION
TOP_N = 100
SAS_THRESHOLD = 4.0
MASTER_KEY = "YourSuperSecretMasterKey"
# Define predefined call limits for each oracle_name. Adjust as needed.
call_limits = {
    "DRD2": float('inf')  # Example: 10 calls for DRD2 oracle. Add other oracles and limits as needed.
}

# App
app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///results.db'
db = SQLAlchemy(app)
socketio = SocketIO(app)

def sum_filter(value):
    return sum(value.values())

app.jinja_env.filters['sum'] = sum_filter