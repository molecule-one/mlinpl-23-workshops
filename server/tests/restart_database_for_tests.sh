#!/bin/bash
# Script to be run before tests
set -e

# rm results.db.tmp || echo 'not found results.db.tmp'
# mv results.db results.db.tmp

# hack: reset both dev and prod default ports
curl -X POST -H "Content-Type: application/json" \
     http://127.0.0.1:5000/reset

curl -X POST -H "Content-Type: application/json" \
     http://127.0.0.1:8000/reset
