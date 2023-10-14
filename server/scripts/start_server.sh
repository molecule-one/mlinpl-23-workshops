#!/bin/bash
# Starts server (for load test or production)
set -e

ulimit -n 1000000

# Pick like 25-50% workers
N_WORKERS=${N_WORKERS:=2}

# Assume 10mins per call max. Gunicorn by default has 60s, which way too low
gunicorn -w ${N_WORKERS} --worker-connections ${N_WORKERS} --worker-class gevent 'server.start:app' --timeout 6000