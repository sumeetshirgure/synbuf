"""Physics check: how much of the decoder work S is explained by the number of faults.

N_err (number of fired DEM error mechanisms) is exactly Poisson-binomial, computable from the
DEM. If S | N_err is concentrated, P(S > x) = sum_n P(N_err = n) P(S > x | N_err = n) can be
evaluated far below the sampled tail, using a parametric model for S | n where data is sparse.
"""
import numpy as np
from scipy import special, stats

from .dem import poisson_binomial_pmf


def dem_error_probs(dem):
    probs = []
    for inst in dem.flattened():
        if inst.type == "error":
            probs.append(inst.args_copy()[0])
    return np.array(probs)


def conditional_stats(S, N, min_count=50):
    """Per value of N: count, mean, median, std of log S, and quantiles."""
    S = np.asarray(S, dtype=np.float64)
    N = np.asarray(N)
    order = np.argsort(N, kind="stable")
    Ns, Ss = N[order], S[order]
    vals, starts, counts = np.unique(Ns, return_index=True, return_counts=True)
    rows = []
    for v, s0, c in zip(vals, starts, counts):
        if c < min_count:
            continue
        x = Ss[s0:s0 + c]
        lx = np.log(np.maximum(x, 1))
        rows.append(dict(n=int(v), count=int(c), mean=float(x.mean()), median=float(np.median(x)),
                         q99=float(np.quantile(x, 0.99)), mean_log=float(lx.mean()), std_log=float(lx.std())))
    return rows


def variance_explained(S, N):
    """Fraction of Var(log S) explained by N (via the conditional mean of log S given N)."""
    lS = np.log(np.maximum(np.asarray(S, dtype=np.float64), 1))
    N = np.asarray(N)
    vals, inv, counts = np.unique(N, return_inverse=True, return_counts=True)
    sums = np.bincount(inv, weights=lS)
    cm = sums / counts
    resid = lS - cm[inv]
    return float(1 - resid.var() / lS.var())


def fit_product_model(S, N, tail_frac=1e-3, min_count=200):
    """Product model S = m(N) * Y with m(n) = exp(a) * n^b and Y independent of N.

    m(n) is fit to the conditional medians; Y = S / m(N) is pooled over all shots with N >= 1 and
    its distribution is represented by the empirical bulk plus a fitted tail (best-AIC model among
    exponential / Weibull / lognormal / Pareto above the (1 - tail_frac) quantile of Y).
    """
    from . import tails

    S = np.asarray(S, dtype=np.float64)
    N = np.asarray(N)
    rows = [r for r in conditional_stats(S, N, min_count=min_count) if r["n"] >= 1]
    n = np.array([r["n"] for r in rows], dtype=np.float64)
    w = np.sqrt(np.array([r["count"] for r in rows], dtype=np.float64))
    med = np.array([r["median"] for r in rows])
    X = np.vstack([np.ones_like(n), np.log(n)]).T
    a, b = np.linalg.lstsq(X * w[:, None], np.log(med) * w, rcond=None)[0]
    m = N >= 1
    Y = S[m] / np.exp(a + b * np.log(N[m]))
    # scale Y so that the tail fit sees values >> 1 (continuity correction is in units of S)
    ta = tails.analyze_tail(Y * 1000.0, tail_frac=tail_frac)
    best = ta["best_aic"]
    ysorted = np.sort(Y)
    # independence check: Y quantiles per n
    indep = []
    for r in rows:
        yn = S[N == r["n"]] / np.exp(a + b * np.log(r["n"]))
        if len(yn) >= 10000:
            indep.append(dict(n=r["n"], q1e_2=float(np.quantile(yn, 1 - 1e-2)),
                              q1e_3=float(np.quantile(yn, 1 - 1e-3)), median=float(np.median(yn))))
    return dict(a=float(a), b=float(b), y_tail=ta, y_best=best, ysorted=ysorted, p_n0=float((N == 0).mean()),
                s_n0_max=float(S[N == 0].max()) if (N == 0).any() else 0.0, independence=indep)


def _y_sf(model, y):
    """P(Y > y): empirical below the tail threshold, fitted model above it."""
    from . import tails

    ta = model["y_tail"]
    u = (ta["u"] + 0.5) / 1000.0
    ys = model["ysorted"]
    y = np.asarray(y, dtype=np.float64)
    emp = (len(ys) - np.searchsorted(ys, y, side="right")) / len(ys)
    best = model["y_best"]
    mod = tails.model_sf(best, ta["fits"][best]["params"], np.maximum(y, u) * 1000.0, u * 1000.0, ta["pu"])
    return np.where(y < u, emp, mod)


def physics_tail(model, pmf, xs):
    """P(S > x) = sum_n P(N_err = n) * P(Y > x / m(n)) (n = 0 contributes only below s_n0_max)."""
    ns = np.arange(1, len(pmf), dtype=np.float64)
    mn = np.exp(model["a"] + model["b"] * np.log(ns))
    out = []
    for x in np.atleast_1d(xs):
        v = float(np.sum(pmf[1:] * _y_sf(model, x / mn)))
        if x < model["s_n0_max"]:
            v += pmf[0]
        out.append(v)
    return np.array(out)


def nerr_pmf(dem, nmax):
    return poisson_binomial_pmf(dem_error_probs(dem), nmax=nmax)
