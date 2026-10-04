#!/usr/bin/env bash
# Secondary run: stop_nconv=5 (collect 5 solutions), reported separately from time to first solution.
cd "$(dirname "$0")/.."
export PYTHONPATH=.
while kill -0 25500 2>/dev/null; do sleep 30; done
python3 -m synbuf.run_relay --codes 72_12_6 --p 0.003 --shots 100000 --stop-nconv 5 --tag _nconv5 >> logs/relay.log 2>&1
python3 -m synbuf.run_relay --codes 144_12_12 --p 0.003 --shots 20000 --stop-nconv 5 --tag _nconv5 >> logs/relay.log 2>&1
