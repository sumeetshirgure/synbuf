#!/usr/bin/env bash
# BB-code data collection with per-point shot budgets (Relay-BP cost is dominated by the heavy tail
# of the iteration count, so the large codes / high-p points get fewer shots). Chunks are resumable.
cd "$(dirname "$0")/.."
export PYTHONPATH=.
run() { python3 -m synbuf.run_relay --codes "$1" --p "$2" --shots "$3" >> logs/relay.log 2>&1; }
# p = 0.003 (reference error rate)
run 90_8_10 0.003 1000000
run 108_8_10 0.003 1000000
run 144_12_12 0.003 500000
run 288_12_18 0.003 200000
# p = 0.002
run 72_12_6 0.002 1000000
run 90_8_10 0.002 1000000
run 108_8_10 0.002 1000000
run 144_12_12 0.002 500000
run 288_12_18 0.002 200000
# p = 0.004
run 72_12_6 0.004 200000
run 90_8_10 0.004 200000
run 108_8_10 0.004 200000
run 144_12_12 0.004 200000
run 288_12_18 0.004 200000
