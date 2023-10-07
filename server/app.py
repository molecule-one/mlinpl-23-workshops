"""Standalone server to keep track of and serve the leaderboard"""
import os
from flask import Flask
from flask_socketio import SocketIO
from flask_sqlalchemy import SQLAlchemy

# CONFIGURATION
TOP_N = 10
SAS_THRESHOLD = 4.0
MASTER_KEY = os.environ.get("MASTER_KEY", "YourSuperSecretMasterKey")
# Define predefined call limits for each oracle_name. Adjust as needed.
call_limits = {
    "DRD2": float('inf')  # Example: 10 calls for DRD2 oracle. Add other oracles and limits as needed.
}
WORKSHOP_ORACLES = ['DRD2', 'JNK3', 'GSK3β']

# App
app = Flask(__name__)
path = os.path.abspath(os.getcwd())
db_path = os.path.join(path, "results.db")
app.config['SQLALCHEMY_DATABASE_URI'] = f'sqlite:///{db_path}'
db = SQLAlchemy(app)
socketio = SocketIO(app)

def sum_filter(value):
    return sum(value.values())

app.jinja_env.filters['sum'] = sum_filter