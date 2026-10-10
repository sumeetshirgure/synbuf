# Third-party patches

synbuf uses instrumented builds of three external decoders. The patches in this directory apply
cleanly to the exact upstream commits pinned below; `scripts/build_decoders.sh` checks out these
commits, applies the patches and installs the packages.

| Package | Upstream repository | Pinned commit | Commit date / subject | Patch |
|---|---|---|---|---|
| PyMatching (Sparse Blossom) | https://github.com/oscarhiggott/PyMatching | `6f63b2b9474ba0fa7e511fe52bffdce858a06984` | 2026-05-21, "Fix CI workflow to trigger on tags (#184)" | `pymatching_op_counters.patch` |
| Relay-BP | https://github.com/trmue/relay | `d185194ba0cb4101ced4340d82b2ee6d42f225f0` | 2026-07-21, "Fix syndrome pruning (#42)" | `relay_bp_stats.patch` |
| ldpc (union-find decoder) | https://github.com/quantumgizmos/ldpc | `d3429964cd4ffe1abfc041c6ec8b8425cb174f40` | 2025-12-08, "Bumping version number after PR #97 (memory leak fix)" | `ldpc_uf_op_counters.patch` |
| pybind11 (PyMatching submodule, unpatched) | https://github.com/pybind/pybind11 | tag `v2.13.6` = `a2e59f0e7065404b44dfe92a28aca47ba1378dc4` | — | none |

The same information is machine-readable in `versions.json`.

## What the patches do

- **pymatching_op_counters.patch** adds per-decode elementary-operation counters to Sparse
  Blossom (`src/pymatching/sparse_blossom/op_counters.h`): priority-queue pops, valid flooder
  events, matcher (MWPM) events, detector-node reschedules and blossom-shatter calls. It also adds a
  `decode_batch_with_op_counts` binding that returns predictions plus an `(N, 5)` counter array.
  Decoding results are unchanged; predictions were checked to be bit-identical to stock
  `decode_batch`.
- **relay_bp_stats.patch** adds a `BPExtraResult::Relay` variant that records the BP iterations up
  to the first converged leg, the number of legs executed and the first converged correction. It
  also adds `ObservableDecoderRunner.decode_stats_batch`, a parallel batch API returning
  `[converged, total_iters, first_conv_iters, legs, obs_best_mask, obs_first_mask]` per shot.
  The decoding algorithm is unchanged.
- **ldpc_uf_op_counters.patch** adds per-decode cluster-operation counters to ldpc's union-find
  decoder (`ldpc::uf::OpCounters` in `src_cpp/union_find.hpp`): singleton clusters created, growth
  rounds, cluster growth steps, bits added, cluster merges, nodes copied by merges, clusters scanned
  for validity, peeling steps and per-cluster GF(2) solves. It also adds
  `UnionFindDecoder.decode_batch_with_op_counts(syndromes, edge_obs)`, which returns predicted
  observable masks plus an `(N, 9)` counter array. The decoding algorithm is unchanged.
  Two upstream quirks matter here. `uf_method` is typed `str`, so the default `False` raises a
  `TypeError`, and any non-empty string (including `"peeling"`) selects matrix-inversion mode. Only
  `uf_method=""` selects peeling. Peeling mode also aborts (`Iterator index out of bounds`) on
  circuit-level surface-code clusters that contain several boundary edges, so synbuf uses
  matrix-inversion mode (`uf_method="inversion"`).
