"""Experiment 1/2 data collection for Relay-BP on bivariate bicycle codes.

Relay-BP runs with an effectively unlimited leg budget (num_sets large) and stop_nconv=1
by default, so S = total BP iterations summed over legs until the first converged leg
("time to first solution"). Per shot we record:
  it_first  - BP iterations up to and including the first converged leg (S)
  it_total  - BP iterations until the stopping criterion (equals it_first when stop_nconv=1)
  legs      - legs executed
  conv      - whether any leg converged (False = censored at the leg limit)
  fail      - logical failure of the reported (best) solution; censored shots count as failures
  fail_first- logical failure of the first converged solution
  ndef/nerr - detection events and fired DEM mechanisms (nerr only on the first nerr_shots shots,
              65535 elsewhere)
Results are written as resumable chunk files.
"""
import argparse
import glob
import json
import os
import time
import zlib

import numpy as np

from .circuits import BB_CODES, bb_memory_circuit
from .relay import COL_CONVERGED, COL_FIRST_ITERS, COL_ITERS, COL_LEGS, COL_OBS_BEST, COL_OBS_FIRST, build_relay, pack_obs

# Relay-BP parameters per code (filled in from scripts/tune_relay.py; defaults = gross-code values
# from Müller et al.).
DEFAULT_PARAMS = dict(gamma0=0.125, gamma_dist_interval=(-0.24, 0.66), pre_iter=80, set_max_iter=60)


def load_params(code, path="results/relay_params.json"):
    params = dict(DEFAULT_PARAMS)
    if os.path.exists(path):
        tuned = json.load(open(path))
        params.update(tuned.get(code, {}))
    params["gamma_dist_interval"] = tuple(params["gamma_dist_interval"])
    return params


def run(code_name, p, shots, outdir, chunk=20_000, num_sets=10_000, stop_nconv=1, seed=0, params=None,
        nerr_shots=200_000):
    code = BB_CODES[code_name]
    params = load_params(code_name) if params is None else params
    os.makedirs(outdir, exist_ok=True)
    meta = dict(code=code_name, p=p, num_sets=num_sets, stop_nconv=stop_nconv, rounds=code.d, **params)
    json.dump(meta, open(os.path.join(outdir, "meta.json"), "w"), indent=1)
    circ = bb_memory_circuit(code, p)
    dem = circ.detector_error_model()
    runner, H, L, pri = build_relay(dem, num_sets=num_sets, stop_nconv=stop_nconv, seed=seed, **params)
    nchunks = (shots + chunk - 1) // chunk
    t0 = time.time()
    for ci in range(nchunks):
        path = os.path.join(outdir, f"part_{ci:05d}.npz")
        if os.path.exists(path):
            continue
        n = min(chunk, shots - ci * chunk)
        sampler = dem.compile_sampler(seed=seed * 1_000_003 + ci + 1)
        want_err = ci * chunk < nerr_shots  # stim's error-returning sampler is slow; subset only
        if want_err:
            det, obs, err = sampler.sample(n, bit_packed=True, return_errors=True)
        else:
            det, obs, _ = sampler.sample(n, bit_packed=True)
        det_u = np.unpackbits(det, axis=1, count=dem.num_detectors, bitorder="little")
        obs_u = np.unpackbits(obs, axis=1, count=dem.num_observables, bitorder="little")
        st = runner.decode_stats_batch(det_u, parallel=True)
        om = pack_obs(obs_u)
        conv = st[:, COL_CONVERGED].astype(bool)
        np.savez(path,
                 it_first=st[:, COL_FIRST_ITERS].astype(np.uint32),
                 it_total=st[:, COL_ITERS].astype(np.uint32),
                 legs=st[:, COL_LEGS].astype(np.uint32),
                 conv=conv,
                 fail=(st[:, COL_OBS_BEST] != om) | ~conv,
                 fail_first=(st[:, COL_OBS_FIRST] != om) | ~conv,
                 ndef=np.bitwise_count(det).sum(axis=1).astype(np.uint16),
                 nerr=(np.bitwise_count(err).sum(axis=1).astype(np.uint16) if want_err
                       else np.full(n, np.iinfo(np.uint16).max, dtype=np.uint16)))
        print(f"{code_name} p={p:g} chunk {ci + 1}/{nchunks} elapsed {time.time() - t0:.0f}s", flush=True)


def load(outdir):
    parts = sorted(glob.glob(os.path.join(outdir, "part_*.npz")))
    if not parts:
        return None
    data = [dict(np.load(f)) for f in parts]
    res = {k: np.concatenate([q[k] for q in data]) for k in data[0]}
    res["meta"] = json.load(open(os.path.join(outdir, "meta.json")))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--codes", nargs="+", required=True)
    ap.add_argument("--p", type=float, nargs="+", required=True)
    ap.add_argument("--shots", type=int, default=1_000_000)
    ap.add_argument("--chunk", type=int, default=20_000)
    ap.add_argument("--num-sets", type=int, default=10_000)
    ap.add_argument("--stop-nconv", type=int, default=1)
    ap.add_argument("--outroot", default="data/relay")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    for p in a.p:
        for code in a.codes:
            outdir = os.path.join(a.outroot, f"bb{code}_p{p:g}{a.tag}")
            run(code, p, a.shots, outdir, chunk=a.chunk, num_sets=a.num_sets, stop_nconv=a.stop_nconv,
                seed=int(p * 1e5) + zlib.crc32(code.encode()) % 1000)
            r = load(outdir)
            print(f"DONE {code} p={p:g} shots={len(r['fail'])} LER={r['fail'].mean():.3e} "
                  f"S mean={r['it_first'].mean():.1f} max={r['it_first'].max()} censored={(~r['conv']).sum()}",
                  flush=True)
