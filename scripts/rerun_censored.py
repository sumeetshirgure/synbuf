"""Re-decode censored Relay-BP shots (no leg converged) with a much larger leg budget.

Regenerates the syndromes deterministically from the chunk seeds used by synbuf.run_relay, picks
up to --max shots that were censored, and decodes them with --num-sets legs. Tells us whether
censored shots eventually converge (power-law tail) or never do (an atom at infinity).
"""
import argparse
import glob
import json
import os
import sys
import time
import zlib

import numpy as np

sys.path.insert(0, ".")
from synbuf.circuits import BB_CODES, bb_memory_circuit  # noqa: E402
from synbuf.relay import COL_CONVERGED, COL_FIRST_ITERS, COL_LEGS, COL_OBS_FIRST, build_relay, pack_obs  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--code", required=True)
ap.add_argument("--p", type=float, required=True)
ap.add_argument("--num-sets", type=int, default=100_000)
ap.add_argument("--max", type=int, default=20)
ap.add_argument("--chunk", type=int, default=20_000)
a = ap.parse_args()

outdir = f"data/relay/bb{a.code}_p{a.p:g}"
meta = json.load(open(os.path.join(outdir, "meta.json")))
seed = int(a.p * 1e5) + zlib.crc32(a.code.encode()) % 1000
code = BB_CODES[a.code]
dem = bb_memory_circuit(code, a.p).detector_error_model()
dets, obss, where = [], [], []
for f in sorted(glob.glob(os.path.join(outdir, "part_*.npz"))):
    ci = int(f[-9:-4])
    conv = np.load(f)["conv"]
    bad = np.nonzero(~conv)[0]
    if not len(bad):
        continue
    n = len(conv)
    det, obs, _ = dem.compile_sampler(seed=seed * 1_000_003 + ci + 1).sample(n, bit_packed=True)
    det_u = np.unpackbits(det, axis=1, count=dem.num_detectors, bitorder="little")
    obs_u = np.unpackbits(obs, axis=1, count=dem.num_observables, bitorder="little")
    for b in bad:
        dets.append(det_u[b]); obss.append(obs_u[b]); where.append((ci, int(b)))
    if len(dets) >= a.max:
        break
dets = np.array(dets[: a.max], dtype=np.uint8)
obss = np.array(obss[: a.max])
params = {k: meta[k] for k in ("gamma0", "gamma_dist_interval", "pre_iter", "set_max_iter")}
runner, *_ = build_relay(dem, num_sets=a.num_sets, stop_nconv=1, seed=777, **params)
t = time.time()
st = runner.decode_stats_batch(dets, parallel=True)
res = dict(code=a.code, p=a.p, num_sets=a.num_sets, original_limit=meta["pre_iter"] + meta["num_sets"] * meta["set_max_iter"],
           n=len(dets), converged=st[:, COL_CONVERGED].tolist(), iters=st[:, COL_FIRST_ITERS].tolist(),
           legs=st[:, COL_LEGS].tolist(), logical_fail=(st[:, COL_OBS_FIRST] != pack_obs(obss)).tolist(),
           where=where[: a.max], seconds=time.time() - t)
print(json.dumps(res))
json.dump(res, open(f"results/relay_censored_rerun_{a.code}_p{a.p:g}.json", "w"), indent=1)
