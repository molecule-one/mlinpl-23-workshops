from server.app import WORKSHOP_ORACLES
from src.server_wrapper import FlaskAppClient

TEST_TOKEN_PREFIX = 'test-' # test-0, test-1, ...

def test_submitting_compounds_to_workshop_oracles():
    """Submits three simple molecules to the server using token test-0, to all workshop oracles."""
    client = FlaskAppClient()
    token = TEST_TOKEN_PREFIX + '0'

    # Example for scoring compounds

    for oracle in WORKSHOP_ORACLES:
        compounds = ["CCCCCCCCC", "CCCCCCCC", "CCCCCC=O"]
        response = client.score_compounds_and_update_leaderboard(compounds, oracle, token)
        print(response)
        assert "metrics" in response
        assert "compound_scores" in response
        assert "compound_sas_scores" in response


def test_random_exploration_gets_reasonable_score():
    """Simple random exploration of ZINC. Should get above >0.5 score on each oracle."""
    pass