#!/usr/bin/env bash
# Build the instrumented decoders used by synbuf (pinned commits: see third_party/README.md).
#  - PyMatching (Sparse Blossom) with elementary-operation counters (decode_batch_with_op_counts)
#  - Relay-BP with first-solution iteration counts and leg statistics (decode_stats_batch)
#  - ldpc union-find decoder with cluster-operation counters (decode_batch_with_op_counts)
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-/tmp/synbuf_build}"
mkdir -p "$WORK" && cd "$WORK"

git clone https://github.com/oscarhiggott/PyMatching.git pm
(cd pm && git checkout 6f63b2b9474ba0fa7e511fe52bffdce858a06984 && git clone --depth 1 -b v2.13.6 https://github.com/pybind/pybind11.git pybind11 \
   && git apply "$HERE/third_party/pymatching_op_counters.patch" && pip install .)

git clone https://github.com/trmue/relay.git relay
(cd relay && git checkout d185194ba0cb4101ced4340d82b2ee6d42f225f0 && git apply "$HERE/third_party/relay_bp_stats.patch" && pip install .)

git clone https://github.com/quantumgizmos/ldpc.git ldpc
(cd ldpc && git checkout d3429964cd4ffe1abfc041c6ec8b8425cb174f40 && git apply "$HERE/third_party/ldpc_uf_op_counters.patch" && pip install .)

pip install stim numpy scipy matplotlib numba
