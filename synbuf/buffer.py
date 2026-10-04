"""Buffer sizing from tail models (Experiment 2 outputs) plus a Lindley-recursion check.

Queue model: one window of syndrome data arrives every `a` work units (deterministic); a server
processes 1 work unit per unit time; window decode work S is i.i.d.; load rho = E[S]/a.
The backlog (workload) W seen by an arriving window obeys W' = max(0, W + S - a). The syndrome
buffer must hold ~W/a + 1 windows; it overflows when W exceeds the buffer's work capacity x.

Leading-order sizing (single big jump, subexponential S; Veraverbeke):
    P(W > x) ~ rho/(1-rho) * Fbar_I(x),  Fbar_I(x) = (1/E[S]) * int_x^inf Fbar(y) dy
D/G/k with m spare servers and out-of-order release (Foss–Korshunov): overflow needs ~m
simultaneous big jumps, P(W > x) ~ C * Fbar_I(x)^m; to leading order x_m(eps) = Fbar_I^{-1}(eps^{1/m}).
The plain quantile Fbar^{-1}(eps) is the "table" leading order used in EXPERIMENTS.md.
"""
import numpy as np
from scipy import integrate, optimize

try:
    import numba
except ImportError:  # pragma: no cover
    numba = None


def integrated_tail_fn(sf, mean, xmax_factor=1e6):
    """Return Fbar_I(x) = (1/mean) * int_x^inf sf(y) dy evaluated numerically in log-space."""
    def FI(x):
        f = lambda t: sf(np.exp(t)) * np.exp(t)
        lo = np.log(x)
        val, _ = integrate.quad(f, lo, lo + np.log(xmax_factor), limit=400, epsabs=0, epsrel=1e-6)
        return val / mean
    return FI


def solve_x(fn, target, x_lo, x_hi=1e30):
    """Solve fn(x) = target for decreasing fn on [x_lo, x_hi] (log-space bisection)."""
    a, b = np.log(x_lo), np.log(x_hi)
    fa = fn(np.exp(a)) - target
    if fa <= 0:
        return x_lo
    for _ in range(200):
        m = 0.5 * (a + b)
        if fn(np.exp(m)) - target > 0:
            a = m
        else:
            b = m
        if b - a < 1e-6:
            break
    return float(np.exp(0.5 * (a + b)))


if numba is not None:
    @numba.njit(cache=True)
    def _lindley(S, a, thresholds, W):
        cnt = np.zeros(len(thresholds), dtype=np.int64)
        for s in S:
            W = W + s - a
            if W < 0:
                W = 0.0
            for i in range(len(thresholds)):
                if W > thresholds[i]:
                    cnt[i] += 1
        return cnt, W
else:  # pragma: no cover
    def _lindley(S, a, thresholds, W):
        cnt = np.zeros(len(thresholds), dtype=np.int64)
        for s in S:
            W = max(0.0, W + s - a)
            cnt += W > thresholds
        return cnt, W


def lindley_overflow(S_samples, rho, thresholds, n_windows, rng=None, block=10_000_000):
    """Simulate D/G/1 by resampling S_samples; return P(W > x) for each threshold x."""
    rng = np.random.default_rng(1) if rng is None else rng
    S_samples = np.asarray(S_samples, dtype=np.float64)
    a = S_samples.mean() / rho
    thresholds = np.asarray(thresholds, dtype=np.float64)
    tot = np.zeros(len(thresholds), dtype=np.int64)
    done = 0
    W = 0.0
    while done < n_windows:
        n = min(block, n_windows - done)
        S = rng.choice(S_samples, size=n, replace=True)
        cnt, W = _lindley(S, a, thresholds, W)
        tot += cnt
        done += n
    return tot / n_windows, a
