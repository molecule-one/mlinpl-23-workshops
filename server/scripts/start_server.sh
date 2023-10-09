#!/bin/bash
# Starts server (for load test or production)
set -e

ulimit -n 1000000

gunicorn -w 4 --worker-connections 4 --worker-class gevent 'server.start:app'