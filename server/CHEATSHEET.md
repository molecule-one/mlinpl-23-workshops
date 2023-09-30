# Useful commands

1. Generate tokens for workshop attendees
```commandline
curl -X POST -H "Content-Type: application/json" \
     -d '{"master_key": "YourSuperSecretMasterKey"}' \
     http://127.0.0.1:5000/generate_tokens
```

## 2 Adding results

CURL:

```commandline
curl -X POST -H "Content-Type: application/json" \
-d '{
    "token": "unique_token_1",
    "metrics": {
        "metric1": 10,
        "metric2": 20
    }
}' \
http://127.0.0.1:5000/add_result

```

You can also use `tests/add_some_results.sh`.