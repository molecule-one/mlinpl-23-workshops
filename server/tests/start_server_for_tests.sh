#!/bin/bash
# Script to be run before tests
set -e

ulimit -n 1000000

PORT=${PORT:=5000}

python server/start.py
