#!/bin/bash
# Script to be run before tests
set -e

conda activate de-novo-workshop
python server/start.py

curl -X POST -H "Content-Type: application/json" \
     http://127.0.0.1:5000/reset
