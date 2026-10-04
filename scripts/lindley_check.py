"""D/G/1 Lindley-recursion check of the leading-order buffer formula, using resampled empirical S.

Compares the simulated P(W > x) with rho/(1-rho) * Fbar_I(x) computed from the same empirical
distribution, at loads rho = 0.5 and 0.8. Writes results/exp2_lindley_<tag>.json.
"""
import json
import sys

import numpy as np

sys.path.insert(0, ".")
from synbuf import buffer as buf  # noqa: E402

path, tag = sys.argv[1], sys.argv[2]
n_windows = int(float(sys.argv[3])) if len(sys.argv) > 3 else 200_000_000
if path.endswith(".npz"):
    S = np.load(path)["S"].astype(np.float64)
else:
    from synbuf.run_relay import load
    S = load(path)["it_first"].astype(np.float64)
med, ES = np.median(S), S.mean()
thr_med = np.array([2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000], dtype=float)
thr = thr_med * med
Ss = np.sort(S)
csum = np.r_[0, np.cumsum(Ss[::-1])][::-1]  # sum of S[i:] for sorted S
out = dict(tag=tag, n_samples=len(S), median=med, mean=ES, thresholds_med=thr_med.tolist(), rho={})
for rho in (0.5, 0.8):
    pr, a = buf.lindley_overflow(S, rho, thr, n_windows=n_windows)
    # empirical integrated tail: E[(S-x)^+]/E[S]
    FI = []
    for x in thr:
        i = np.searchsorted(Ss, x, side="right")
        FI.append((csum[i] - x * (len(Ss) - i)) / len(Ss) / ES)
    out["rho"][str(rho)] = dict(a=a, P_overflow_sim=pr.tolist(),
                                veraverbeke=[rho / (1 - rho) * f for f in FI], n_windows=n_windows)
    print(tag, rho, [f"{t:g}:{p:.2e}/{v:.2e}" for t, p, v in zip(thr_med, pr, out["rho"][str(rho)]["veraverbeke"])], flush=True)
json.dump(out, open(f"results/exp2_lindley_{tag}.json", "w"), indent=1)
