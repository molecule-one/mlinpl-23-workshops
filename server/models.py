"""Standalone server to keep track of and serve the leaderboard"""
from datetime import datetime

from server.app import db


class User(db.Model):
    id = db.Column(db.String(80), primary_key=True)
    oracle_calls = db.Column(db.PickleType, nullable=False)  # Dictionary with {oracle_name: count}


class Result(db.Model):
    id = db.Column(db.String(80), primary_key=True)
    metrics = db.Column(db.PickleType, nullable=False)


class Token(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    token = db.Column(db.String(128), unique=True, nullable=False)
    creation_date = db.Column(db.DateTime, default=datetime.utcnow)

    @classmethod
    def check_valid_token(cls, token_str):
        return cls.query.filter_by(token=token_str).first() is not None
