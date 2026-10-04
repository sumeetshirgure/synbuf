#!/usr/bin/env bash
# Full data-collection pipeline: MWPM sweep, then Relay-BP tuning and BB data collection.
cd "$(dirname "$0")/.."
export PYTHONPATH=.
python3 -m synbuf.run_mwpm --d 3 5 7 9 11 13 15 --p 0.001 0.002 0.003 0.005 --shots 10000000 >> logs/mwpm.log 2>&1
python3 scripts/tune_relay.py > logs/tune.log 2>&1
for p in 0.003 0.002 0.004; do
  python3 -m synbuf.run_relay --codes 72_12_6 90_8_10 108_8_10 144_12_12 288_12_18 --p $p --shots 1000000 >> logs/relay.log 2>&1
done
