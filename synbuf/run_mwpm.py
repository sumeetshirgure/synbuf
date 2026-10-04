"""Experiment 1/2 data collection for MWPM (instrumented PyMatching) on the rotated surface code.

For every shot records the flooder op counters and the uncapped logical failure indicator.
Decodes always run to completion.

Shots are split into batches of `chunk` shots that are decoded by a pool of `workers` processes.
Batch i is sampled with a seed derived from (seed, i), so results do not depend on `workers`.
"""
import argparse
import multiprocessing as mp
import os
import time

import numpy as np
import pymatching

from .circuits import surface_code_circuit


MAX_SHOTS = 1_000_000

_dem = None
_graph = None


def _init(d, p):
    global _dem, _graph
    _dem = surface_code_circuit(d, p).detector_error_model(decompose_errors=True)
    _graph = pymatching.Matching.from_detector_error_model(_dem)._matching_graph
    if not hasattr(_graph, "decode_batch_with_op_counts"):
        raise RuntimeError("pymatching is not the synbuf-instrumented build")


def _batch(args):
    n, seed = args
    det, obs, _ = _dem.compile_sampler(seed=seed).sample(n, bit_packed=True)
    pred, cnt = _graph.decode_batch_with_op_counts(det)
    return dict(
        pq=cnt[:, 0].astype(np.uint32),
        flood=cnt[:, 1].astype(np.uint32),
        mwpm=cnt[:, 2].astype(np.uint32),
        resched=cnt[:, 3].astype(np.uint32),
        fail=((pred[:, 0] ^ obs[:, 0]) & 1).astype(bool),
    )


def run(d, p, shots, out_path, workers=4, chunk=50_000, seed=0):
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
    res["S"] = (res["pq"].astype(np.uint64) + res["resched"]).astype(np.uint32)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    np.savez(out_path, d=d, p=p, shots=shots, seconds=time.time() - t0, **res)
    return res, time.time() - t0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=int, nargs="+", required=True)
    ap.add_argument("--p", type=float, nargs="+", required=True)
    ap.add_argument("--shots", type=int, default=MAX_SHOTS)
    ap.add_argument("--chunk", type=int, default=50_000)
    ap.add_argument("--workers", "-N", type=int, default=4, help="number of worker processes")
    ap.add_argument("--outdir", default="data/mwpm")
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
