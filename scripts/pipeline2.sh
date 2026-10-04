#!/usr/bin/env bash
# Reduced-grid tuning for the two largest codes, then BB data collection.
cd "$(dirname "$0")/.."
export PYTHONPATH=.
python3 scripts/tune_relay.py 144_12_12 --gamma0 0.2 0.3 --interval -0.12 0.5 >> logs/tune.log 2>&1
python3 scripts/tune_relay.py 288_12_18 --gamma0 0.125 0.2 0.3 --interval -0.12 0.5 --p 0.005 --shots 3000 >> logs/tune.log 2>&1
for p in 0.003 0.002 0.004; do
  python3 -m synbuf.run_relay --codes 72_12_6 90_8_10 108_8_10 144_12_12 288_12_18 --p $p --shots 1000000 >> logs/relay.log 2>&1
done
