"""Experiment 1/2 data collection for MWPM (instrumented PyMatching) on the rotated surface code.

For every shot records: flooder op counters, N_def (detection events) and the uncapped logical
failure indicator; N_err (fired DEM mechanisms) is recorded on the first `nerr_shots` shots
(NERR_MISSING elsewhere). Decodes always run to completion.
"""
import argparse
import multiprocessing as mp
import os
import time

import numpy as np
import pymatching

from .circuits import surface_code_circuit


NERR_MISSING = np.iinfo(np.uint16).max


def _worker(args):
    d, p, shots, chunk, seed, nerr_shots = args
    circ = surface_code_circuit(d, p)
    dem = circ.detector_error_model(decompose_errors=True)
    m = pymatching.Matching.from_detector_error_model(dem)
    g = m._matching_graph
    if not hasattr(g, "decode_batch_with_op_counts"):
        raise RuntimeError("pymatching is not the synbuf-instrumented build")
    sampler = dem.compile_sampler(seed=seed)
    out = {k: [] for k in ("pq", "resched", "mwpm", "flood", "ndef", "nerr", "fail")}
    done = 0
    while done < shots:
        n = min(chunk, shots - done)
        # Stim's error-returning sampler is ~30x slower, so N_err is recorded on a subset only.
        want_err = done < nerr_shots
        if want_err:
            det, obs, err = sampler.sample(n, bit_packed=True, return_errors=True)
        else:
            det, obs, _ = sampler.sample(n, bit_packed=True)
        pred, cnt = g.decode_batch_with_op_counts(det)
        out["pq"].append(cnt[:, 0].astype(np.uint32))
        out["flood"].append(cnt[:, 1].astype(np.uint32))
        out["mwpm"].append(cnt[:, 2].astype(np.uint32))
        out["resched"].append(cnt[:, 3].astype(np.uint32))
        out["ndef"].append(np.bitwise_count(det).sum(axis=1).astype(np.uint16))
        if want_err:
            out["nerr"].append(np.bitwise_count(err).sum(axis=1).astype(np.uint16))
        else:
            out["nerr"].append(np.full(n, NERR_MISSING, dtype=np.uint16))
        out["fail"].append(((pred[:, 0] ^ obs[:, 0]) & 1).astype(bool))
        done += n
    return {k: np.concatenate(v) for k, v in out.items()}


def run(d, p, shots, out_path, workers=4, chunk=50_000, seed=0, nerr_shots=1_000_000):
    per = [shots // workers + (1 if i < shots % workers else 0) for i in range(workers)]
    tasks = [(d, p, per[i], chunk, seed * 1000 + i + 1, nerr_shots // workers) for i in range(workers)]
    t0 = time.time()
    with mp.Pool(workers) as pool:
        parts = pool.map(_worker, tasks)
    res = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    res["S"] = (res["pq"].astype(np.uint64) + res["resched"]).astype(np.uint32)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    np.savez(out_path, d=d, p=p, shots=shots, seconds=time.time() - t0, **res)
    return res, time.time() - t0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=int, nargs="+", required=True)
    ap.add_argument("--p", type=float, nargs="+", required=True)
    ap.add_argument("--shots", type=int, default=10_000_000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--outdir", default="data/mwpm")
    a = ap.parse_args()
    for p in a.p:
        for d in a.d:
            path = os.path.join(a.outdir, f"surface_d{d}_p{p:g}.npz")
            if os.path.exists(path):
                print("skip", path, flush=True)
                continue
            res, dt = run(d, p, a.shots, path, workers=a.workers, seed=int(d * 100 + p * 1e5))
            print(f"d={d} p={p:g} shots={a.shots} {dt:.0f}s LER={res['fail'].mean():.3e} "
                  f"S mean={res['S'].mean():.1f} max={res['S'].max()} ndef mean={res['ndef'].mean():.1f}",
                  flush=True)
