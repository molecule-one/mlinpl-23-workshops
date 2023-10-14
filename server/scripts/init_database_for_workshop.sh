#!/bin/bash
# Initializes database for workshop and prints tokens for participants
set -e

PORT=${PORT:=5000}

# restarts database and adds test tokens
curl -X POST -H "Content-Type: application/json" \
     http://127.0.0.1:${PORT}/reset

# add tokens for participants and print them
curl -X POST -H "Content-Type: application/json" \
     -d '{"master_key": "YourSuperSecretMasterKey"}' \
     http://127.0.0.1:${PORT}/generate_tokens