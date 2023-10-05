from server.app import WORKSHOP_ORACLES
from src.server_wrapper import FlaskAppClient

TEST_TOKEN_PREFIX = 'test-' # test-0, test-1, ...

def test_submitting_compounds_to_workshop_oracles():
    client = FlaskAppClient()
    token = "test-0"

    # Example for scoring compounds

    for oracle in WORKSHOP_ORACLES:
        compounds = ["CC", "CCC", "CCC=O"]
        response = client.score_compounds_and_update_leaderboard(compounds, oracle, token)
        print(response)
        assert "metrics" in response
        assert "compound_scores" in response
        assert "compound_sas_scores" in response
