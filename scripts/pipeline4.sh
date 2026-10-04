#!/usr/bin/env bash
# Main BB collection (resumable), then the stop_nconv=5 secondary run.
cd "$(dirname "$0")/.."
export PYTHONPATH=.
scripts/pipeline3.sh
python3 -m synbuf.run_relay --codes 72_12_6 --p 0.003 --shots 100000 --stop-nconv 5 --tag _nconv5 >> logs/relay.log 2>&1
python3 -m synbuf.run_relay --codes 144_12_12 --p 0.003 --shots 20000 --stop-nconv 5 --tag _nconv5 >> logs/relay.log 2>&1
