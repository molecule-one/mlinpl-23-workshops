"""Simple wrapper for app"""
import json
from rich.console import Console
from typing import List
import requests

from src.utils import Oracles

APP_URL = "http://localhost:5000"

class FlaskAppClient:
    ERROR_KEY = "error"
    TRACEBACK_KEY = "traceback"

    def __init__(self, base_url="http://localhost:5000"):
        self.base_url = base_url
        self.console = Console()

    def _handle_response(self, response):
        try:
            response_data = response.json()
        except json.JSONDecodeError:
            self.console.print("[red]Failed to parse server response as JSON[/red]")
            response.raise_for_status()  # This will raise an HTTPError if the HTTP request returned an unsuccessful status code.

        if response.status_code == 200:
            return response_data
        else:
            error = response_data.get(self.ERROR_KEY, 'Unknown error')
            tb = response_data.get(self.TRACEBACK_KEY, None)
            self.console.print(f"[red]Server error: {error}[/red]")
            if tb:
                self.console.print(f"[yellow]{tb}[/yellow]")
            response.raise_for_status()

    def add_result(self, user_token, metrics):
        payload = {"token": user_token, "metrics": metrics}
        response = requests.post(f"{self.base_url}/add_result", json=payload)
        return self._handle_response(response)

    def evaluate_and_add_result(self, compounds, oracle_name, user_token):
        payload = {
            "compounds": ",".join(compounds),
            "oracle_name": oracle_name,
            "token": user_token
        }
        response = requests.post(f"{self.base_url}/evaluate_and_add_result", json=payload)
        return self._handle_response(response)

    def score_compound(self, compounds, oracle_name, user_token):
        payload = {
            "compound": compounds,
            "oracle_name": oracle_name,
            "token": user_token
        }
        response = requests.post(f"{self.base_url}/score_compound", json=payload)
        return self._handle_response(response)

    def score_compounds(self, compounds: List, oracle_name: Oracles, token):
        return [self.score_compound(c, oracle_name, token)['score'] for c in compounds]


# Usage Example:
if __name__ == "__main__":
    client = FlaskAppClient()

    # Example for adding results
    token = "your_unique_token"
    metrics = {"metric1": 85, "metric2": 90}
    response = client.add_result(token, metrics)
    print(response)

    # Example for scoring compounds
    compounds = ["CC", "CCC"]
    oracle_name = "DRD2"
    response = client.score_compounds(compounds, oracle_name, token)
    print(response)

    # Example of error handling
    compounds = ["Cxxxxx"]
    oracle_name = "DRD2"
    response = client.score_compounds(compounds, oracle_name, token)
    print(response)