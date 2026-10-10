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

All counters in the PyMatching and ldpc patches count operations of bounded cost, i.e. O(1), or
bounded by the (fixed) maximum degree of the decoding graph. Under that assumption a decode's running
time is within constant factors of the sum of its counters, and the tail class of the sum equals that
of its heaviest-tailed component, whatever the per-counter costs.

- **pymatching_op_counters.patch** adds per-decode elementary-operation counters to Sparse
  Blossom (`src/pymatching/sparse_blossom/op_counters.h`) and a `decode_batch_with_op_counts`
  binding that returns predictions plus an `(N, 11)` counter array:
  - dual (flooder): priority-queue pops (including stale events), radix-heap bucket-refill moves,
    valid flood events, detector-node reschedules and the neighbour scans inside them;
  - primal (matcher): matcher events, alternating-tree operations (tree rotations, common-ancestor
    walks, tree shattering into matches, path pruning, child-list edits), region operations
    (blossom-cycle construction and shattering, recursion over blossom children, blossom-parent
    chain walks) and detector-node visits by region-wide updates (growth-rate changes, blossom
    wrapping/unwrapping, shell cleanup);
  - match extraction: top-level and recursive blossom-shatter steps.

  Decoding results are unchanged; predictions were checked to be bit-identical to stock
  `decode_batch`. `pq_pops` and `node_reschedules` are counted exactly as in the earlier version of
  this patch, so the earlier metric (their sum, stored as `S_legacy`) is reproducible.
- **relay_bp_stats.patch** adds a `BPExtraResult::Relay` variant that records the BP iterations up
  to the first converged leg, the number of legs executed and the first converged correction. It
  also adds `ObservableDecoderRunner.decode_stats_batch`, a parallel batch API returning
  `[converged, total_iters, first_conv_iters, legs, obs_best_mask, obs_first_mask]` per shot.
  The decoding algorithm is unchanged.
- **ldpc_uf_op_counters.patch** instruments ldpc's union-find decoder for graphs (peeling mode) and
  fixes its peeler:
  - Counters (`ldpc::uf::OpCounters` in `src_cpp/union_find.hpp`): clusters created, growth rounds,
    cluster growth calls, check-row entries scanned for candidate bits, bit additions and their
    column visits, merges and the nodes they copy, per-round validity checks (O(1) parity and
    boundary flag per cluster), size-sort comparisons, spanning-tree construction visits and
    peeling steps. Matrix mode only counts its per-cluster GF(2) solves and is not used by synbuf.
  - `UnionFindDecoder.decode_batch_with_op_counts(syndromes, edge_obs)` returns predicted
    observable masks plus an `(N, 13)` counter array.
  - Peeling fix: upstream's spanning-tree peeler treats the virtual boundary node as a leaf even
    when several tree bits attach to it, may peel a boundary bit that is not in the tree, and reads
    an uninitialized bit index when a leaf has no tree bit left. On circuit-level surface-code
    clusters it aborts with `Iterator index out of bounds`. The patch replaces it with the standard
    Delfosse–Zémor peeler on a BFS spanning tree, rooted at the virtual boundary node when the
    cluster touches the boundary. Cluster growth, merging and validity are unchanged. Corrections
    were checked to reproduce the syndrome exactly.

  Upstream quirk: `uf_method` is typed `str`, so the default `False` raises a `TypeError`, and any
  non-empty string (including `"peeling"`) selects matrix-inversion mode. Only `uf_method=""`
  selects peeling.
