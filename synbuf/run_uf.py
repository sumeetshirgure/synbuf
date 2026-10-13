"""Union-find decoder data collection (instrumented ldpc UnionFindDecoder) on the rotated surface code.

For every shot records the union-find cluster-operation counters and the uncapped logical failure
indicator. Decodes always run to completion.

Decoder: ldpc's UnionFindDecoder (Roffe et al.) in graph (peeling) mode, unweighted growth (all
neighbouring bits per step). Validity is the O(1) parity + boundary test; the correction comes from
peeling a BFS spanning tree of each cluster. synbuf's ldpc patch replaces upstream's spanning-tree
peeler, which aborts on circuit-level surface-code clusters (see third_party/README.md).

The decoding graph is the one MWPM uses: the DEM is decomposed into graphlike components
(``decompose_errors=True``); each component is a column of H (weight 2, or weight 1 for a boundary
edge). Parallel edges with the same detectors and observables are merged by XOR-combining their
probabilities; parallel edges that differ in observables keep the more likely one.

Work metric: S = sum of all counters in COUNTERS. Every counter counts bounded-cost operations
(hash-set operations or sparse-matrix entry visits): cluster creation, growth (row scans, bit
additions, column visits), merges (nodes copied), validity checks, sorting, spanning-tree
construction and peeling.

Shots are split into batches of `chunk` shots that are decoded by a pool of `workers` processes.
Batch i is sampled with a seed derived from (seed, i), so results do not depend on `workers`.
"""
import argparse
import multiprocessing as mp
import os
import time

import numpy as np
import scipy.sparse as sp

from .circuits import surface_code_circuit


MAX_SHOTS = 1_000_000
COUNTERS = ["clusters_created", "growth_rounds", "cluster_grows", "growth_row_visits", "bit_add_calls",
            "growth_col_visits", "merges", "merge_moves", "cluster_scans", "sort_compares", "peel_visits",
            "peel_steps", "cluster_solves"]  # cluster_solves: matrix mode only (always 0 here)

_dem = None
_decoder = None
_edge_obs = None


def dem_to_graph(dem):
    """Return (H, edge_obs, probs) for a decomposed (graphlike) DEM.

    H: detectors x edges (csr, uint8); edge_obs[j]: observable bitmask flipped by edge j.
    """
    merged = {}  # dets -> (p, obs)
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        p = inst.args_copy()[0]
        comps, dets, obs = [], [], 0
        for t in inst.targets_copy():
            if t.is_separator():
                comps.append((dets, obs))
                dets, obs = [], 0
            elif t.is_relative_detector_id():
                dets.append(t.val)
            elif t.is_logical_observable_id():
                obs ^= 1 << t.val
        comps.append((dets, obs))
        for dets, obs in comps:
            key = tuple(sorted(dets))
            if not key:
                continue
            if len(key) > 2:
                raise ValueError(f"non-graphlike component {key}")
            if key in merged:
                q, o = merged[key]
                if o == obs:
                    merged[key] = (q * (1 - p) + p * (1 - q), o)
                elif p > q:
                    merged[key] = (p, obs)
            else:
                merged[key] = (p, obs)
    keys = list(merged)
    rows = [v for k in keys for v in k]
    cols = [j for j, k in enumerate(keys) for _ in k]
    H = sp.csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(dem.num_detectors, len(keys)))
    edge_obs = np.array([merged[k][1] for k in keys], dtype=np.uint64)
    probs = np.array([merged[k][0] for k in keys])
    return H, edge_obs, probs


def make_decoder(dem):
    from ldpc.union_find_decoder import UnionFindDecoder

    H, edge_obs, _ = dem_to_graph(dem)
    dec = UnionFindDecoder(H, uf_method="")  # "" selects peeling (any non-empty string selects inversion)
    if not hasattr(dec, "decode_batch_with_op_counts"):
        raise RuntimeError("ldpc is not the synbuf-instrumented build")
    return dec, edge_obs


def _init(d, p):
    global _dem, _decoder, _edge_obs
    _dem = surface_code_circuit(d, p).detector_error_model(decompose_errors=True)
    _decoder, _edge_obs = make_decoder(_dem)


def _batch(args):
    n, seed = args
    det, obs, _ = _dem.compile_sampler(seed=seed).sample(n)
    pred, cnt = _decoder.decode_batch_with_op_counts(det, _edge_obs)
    if cnt.shape[1] != len(COUNTERS):
        raise RuntimeError("ldpc build does not match synbuf's op-counter patch")
    res = {k: cnt[:, i].astype(np.uint32) for i, k in enumerate(COUNTERS)}
    res["fail"] = (pred & 1).astype(bool) ^ obs[:, 0]
    return res


def run(d, p, shots, out_path, workers=4, chunk=10_000, seed=0):
    nchunks = (shots + chunk - 1) // chunk
    tasks = [(min(chunk, shots - ci * chunk), seed * 1_000_003 + ci + 1) for ci in range(nchunks)]
    t0 = time.time()
    if workers <= 1:
        _init(d, p)
        parts = [_batch(t) for t in tasks]
    else:
        with mp.Pool(min(workers, nchunks), initializer=_init, initargs=(d, p)) as pool:
            parts = pool.map(_batch, tasks, chunksize=1)
    res = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    res["S"] = sum(res[k].astype(np.uint64) for k in COUNTERS).astype(np.uint32)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    np.savez(out_path, d=d, p=p, shots=shots, seconds=time.time() - t0, **res)
    return res, time.time() - t0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=int, nargs="+", required=True)
    ap.add_argument("--p", type=float, nargs="+", required=True)
    ap.add_argument("--shots", type=int, default=MAX_SHOTS)
    ap.add_argument("--chunk", type=int, default=10_000)
    ap.add_argument("--workers", "-N", type=int, default=4, help="number of worker processes")
    ap.add_argument("--outdir", default="data/uf")
    a = ap.parse_args()
    if a.shots > MAX_SHOTS:
        ap.error(f"--shots is capped at {MAX_SHOTS}")
    for p in a.p:
        for d in a.d:
            path = os.path.join(a.outdir, f"surface_d{d}_p{p:g}.npz")
            if os.path.exists(path):
                print("skip", path, flush=True)
                continue
            res, dt = run(d, p, a.shots, path, workers=a.workers, chunk=a.chunk, seed=int(d * 100 + p * 1e5))
            print(f"d={d} p={p:g} shots={a.shots} {dt:.0f}s LER={res['fail'].mean():.3e} "
                  f"S mean={res['S'].mean():.1f} max={res['S'].max()}", flush=True)
