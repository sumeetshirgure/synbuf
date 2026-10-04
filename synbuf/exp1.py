"""Experiment 1: cap vs LER floor, computed offline from uncapped per-shot (S, F_inf) data.

P_L(c)      = P(F_inf, S <= c) + P(S > c)
dP_L(c)     = P(S > c) * [1 - P(F_inf | S > c)]
"""
import numpy as np
from scipy import stats


def wilson(k, n, z=1.0):
    """Wilson score interval (default 1 sigma)."""
    if n == 0:
        return (np.nan, np.nan)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (max(c - h, 0.0), min(c + h, 1.0))


def cap_stats(S, fail, c):
    S = np.asarray(S)
    fail = np.asarray(fail, dtype=bool)
    n = len(S)
    slow = S > c
    k_slow = int(slow.sum())
    k_fail_fast = int((fail & ~slow).sum())
    k_fail_slow = int((fail & slow).sum())
    k_PL = k_fail_fast + k_slow
    return dict(
        c=float(c), n=n,
        p_slow=k_slow / n, p_slow_ci=wilson(k_slow, n),
        PL=k_PL / n, PL_ci=wilson(k_PL, n),
        pL0=int(fail.sum()) / n, pL0_ci=wilson(int(fail.sum()), n),
        dPL=(k_slow - k_fail_slow) / n,
        p_fail_given_slow=(k_fail_slow / k_slow) if k_slow else np.nan,
        p_fail_given_slow_ci=wilson(k_fail_slow, k_slow) if k_slow else (np.nan, np.nan),
        k_slow=k_slow, k_fail_slow=k_fail_slow, k_fail=int(fail.sum()),
    )


def cap_curve(S, fail, cs):
    """Vectorized P_L(c), P(S>c), P(F|S>c) over an array of caps."""
    S = np.asarray(S)
    fail = np.asarray(fail, dtype=bool)
    n = len(S)
    order = np.argsort(S, kind="stable")
    Ss = S[order]
    fs = fail[order]
    cum_fail = np.r_[0, np.cumsum(fs)]  # failures among the first i sorted shots
    idx = np.searchsorted(Ss, cs, side="right")  # number with S <= c
    k_slow = n - idx
    k_fail_fast = cum_fail[idx]
    k_fail_slow = cum_fail[-1] - k_fail_fast
    with np.errstate(invalid="ignore", divide="ignore"):
        pfs = np.where(k_slow > 0, k_fail_slow / np.maximum(k_slow, 1), np.nan)
    return dict(c=np.asarray(cs, dtype=float), p_slow=k_slow / n, PL=(k_fail_fast + k_slow) / n,
                pL0=cum_fail[-1] / n, p_fail_given_slow=pfs, k_slow=k_slow)
