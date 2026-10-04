# Third-party patches

synbuf uses instrumented builds of two external decoders. The patches in this directory apply
cleanly to the exact upstream commits pinned below; `scripts/build_decoders.sh` checks out these
commits, applies the patches and installs the packages.

| Package | Upstream repository | Pinned commit | Commit date / subject | Patch |
|---|---|---|---|---|
| PyMatching (Sparse Blossom) | https://github.com/oscarhiggott/PyMatching | `6f63b2b9474ba0fa7e511fe52bffdce858a06984` | 2026-05-21, "Fix CI workflow to trigger on tags (#184)" | `pymatching_op_counters.patch` |
| Relay-BP | https://github.com/trmue/relay | `d185194ba0cb4101ced4340d82b2ee6d42f225f0` | 2026-07-21, "Fix syndrome pruning (#42)" | `relay_bp_stats.patch` |
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
