"""Light per-code retuning of Relay-BP memory strengths (gamma0, gamma distribution interval).

Relay-BP's published parameters were tuned for the gross code [[144,12,12]]. For each BB code we
evaluate a small grid at a physical error rate where the code's LER is measurable, using the
deployed-style leg budget (300 legs, stop_nconv=1), and keep the configuration with the lowest LER
(ties broken by mean iterations). Writes results/relay_params.json and results/relay_tuning.csv.
"""
import csv
import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from synbuf.circuits import BB_CODES, bb_memory_circuit  # noqa: E402
from synbuf.relay import COL_CONVERGED, COL_FIRST_ITERS, COL_OBS_FIRST, build_relay, pack_obs  # noqa: E402

TUNE = {  # code: (p, shots)
    "72_12_6": (0.004, 20000),
    "90_8_10": (0.004, 12000),
    "108_8_10": (0.004, 12000),
    "144_12_12": (0.005, 8000),
    "288_12_18": (0.006, 3000),
}
GAMMA0 = [0.05, 0.125, 0.2, 0.3]
INTERVALS = [(-0.24, 0.66), (-0.12, 0.5), (-0.36, 0.8)]

codes = sys.argv[1:] or list(TUNE)
os.makedirs("results", exist_ok=True)
rows = []
best = json.load(open("results/relay_params.json")) if os.path.exists("results/relay_params.json") else {}
for code in codes:
    p, N = TUNE[code]
    dem = bb_memory_circuit(BB_CODES[code], p).detector_error_model()
    det, obs, _ = dem.compile_sampler(seed=12345).sample(N)
    det = det.astype(np.uint8)
    om = pack_obs(obs)
    results = []
    for g0, iv in itertools.product(GAMMA0, INTERVALS):
        runner, *_ = build_relay(dem, gamma0=g0, gamma_dist_interval=iv, num_sets=300, stop_nconv=1)
        t = time.time()
        st = runner.decode_stats_batch(det, parallel=True)
        fail = (st[:, COL_OBS_FIRST] != om) | (st[:, COL_CONVERGED] == 0)
        r = dict(code=code, p=p, shots=N, gamma0=g0, lo=iv[0], hi=iv[1], failures=int(fail.sum()),
                 ler=float(fail.mean()), mean_iter=float(st[:, COL_FIRST_ITERS].mean()),
                 censored=int((st[:, COL_CONVERGED] == 0).sum()), seconds=time.time() - t)
        print(r, flush=True)
        results.append(r)
        rows.append(r)
    b = min(results, key=lambda r: (r["failures"], r["mean_iter"]))
    best[code] = dict(gamma0=b["gamma0"], gamma_dist_interval=[b["lo"], b["hi"]])
    print("BEST", code, best[code], flush=True)
    json.dump(best, open("results/relay_params.json", "w"), indent=1)

write_header = not os.path.exists("results/relay_tuning.csv")
with open("results/relay_tuning.csv", "a", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    if write_header:
        w.writeheader()
    w.writerows(rows)
