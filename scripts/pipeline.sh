#!/usr/bin/env bash
# Data-collection pipeline: MWPM sweep, Relay-BP tuning (only if results/relay_params.json is
# missing), BB-code collection, the stop_nconv=5 secondary run and the union-find (ldpc) surface-code
# sweep. All runs are resumable.
#
# Usage: scripts/pipeline.sh N [stage ...]
#   N      number of worker processes per run
#   stage  any of: mwpm tune relay nconv5 uf (default: all, in that order)
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=.
N="${1:?usage: scripts/pipeline.sh N [mwpm|tune|relay|nconv5|uf ...]}"
shift
STAGES="${*:-mwpm tune relay nconv5 uf}"
PS="0.0005 0.001 0.003"
SHOTS=1000000
mkdir -p logs

for stage in $STAGES; do
  case "$stage" in
    mwpm)
      python3 -m synbuf.run_mwpm --d 3 7 11 15 --p $PS --shots $SHOTS --workers "$N" >> logs/mwpm.log 2>&1 ;;
    tune)
      [ -f results/relay_params.json ] || python3 scripts/tune_relay.py > logs/tune.log 2>&1 ;;
    relay)
      python3 -m synbuf.run_relay --codes 72_12_6 90_8_10 108_8_10 144_12_12 --p $PS \
        --shots $SHOTS --workers "$N" >> logs/relay.log 2>&1 ;;
    nconv5)
      python3 -m synbuf.run_relay --codes 72_12_6 --p 0.003 --shots 100000 --stop-nconv 5 --tag _nconv5 \
        --workers "$N" >> logs/relay.log 2>&1
      python3 -m synbuf.run_relay --codes 144_12_12 --p 0.003 --shots 20000 --stop-nconv 5 --tag _nconv5 \
        --workers "$N" >> logs/relay.log 2>&1 ;;
    uf)
      python3 -m synbuf.run_uf --d 11 15 --p 0.001 0.003 --shots $SHOTS --workers "$N" >> logs/uf.log 2>&1 ;;
    *)
      echo "unknown stage: $stage" >&2; exit 1 ;;
  esac
done
