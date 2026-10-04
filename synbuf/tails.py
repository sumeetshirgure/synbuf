"""Tail-class characterization of decoder work distributions (Experiment 2).

All fits are conditional on S > u (peaks over threshold) and use left-truncated likelihoods,
so the bulk of the distribution does not influence the tail class. Work counts are integers;
we use the continuity-corrected value S + U(0,1) jitter-free midpoint (S + 0.5) for continuous
models, which is immaterial when u is in the hundreds or more.
"""
import numpy as np
from scipy import optimize, special, stats


# ------------------------------------------------------------------ empirical helpers
def survival(S):
    """Return (x, P(S > x)) at the unique values of S (x sorted ascending)."""
    S = np.sort(np.asarray(S))
    n = len(S)
    x, idx = np.unique(S, return_index=True)
    # number of samples strictly greater than x = n - (index of last occurrence + 1)
    last = np.r_[idx[1:], n] - 1
    return x, (n - last - 1) / n


def survival_at(S, c):
    S = np.sort(np.asarray(S))
    return (len(S) - np.searchsorted(S, c, side="right")) / len(S)


def quantile_tail(S, q):
    """Smallest c with P(S > c) <= q (empirical)."""
    S = np.sort(np.asarray(S))
    n = len(S)
    k = int(np.floor(q * n))  # allowed number of exceedances
    if k >= n:
        return S[0] - 1
    return S[n - k - 1]


def mean_excess(S, us):
    S = np.asarray(S, dtype=np.float64)
    out = []
    for u in us:
        e = S[S > u]
        out.append((e - u).mean() if len(e) else np.nan)
    return np.array(out)


def hill(S, ks):
    """Hill estimator of the tail index alpha using the top-k order statistics."""
    Xs = np.sort(np.asarray(S, dtype=np.float64))[::-1]
    out = []
    for k in ks:
        if k + 1 >= len(Xs) or Xs[k] <= 0:
            out.append(np.nan)
            continue
        out.append(1.0 / np.mean(np.log(Xs[:k]) - np.log(Xs[k])))
    return np.array(out)


# ------------------------------------------------------------------ truncated models
# Each model: logpdf(x, *theta) and logsf(x, *theta) of the *untruncated* distribution,
# parameter transforms to an unconstrained space, and a starting point.
def _exp_logpdf(x, th):
    (lt,) = th
    t = np.exp(lt)
    return lt - t * x


def _exp_logsf(x, th):
    (lt,) = th
    return -np.exp(lt) * x


def _weib_logpdf(x, th):
    ll, lb = th
    lam, b = np.exp(ll), np.exp(lb)
    z = x / lam
    return lb - ll + (b - 1) * np.log(z) - z ** b


def _weib_logsf(x, th):
    ll, lb = th
    return -(x / np.exp(ll)) ** np.exp(lb)


def _ln_logpdf(x, th):
    mu, ls = th
    s = np.exp(ls)
    z = (np.log(x) - mu) / s
    return -0.5 * z * z - np.log(x) - ls - 0.5 * np.log(2 * np.pi)


def _ln_logsf(x, th):
    mu, ls = th
    return special.log_ndtr(-(np.log(x) - mu) / np.exp(ls))


def _par_logpdf(x, th):
    (la,) = th
    a = np.exp(la)
    return la - (a + 1) * np.log(x)  # x0 cancels in truncated likelihood


def _par_logsf(x, th):
    (la,) = th
    return -np.exp(la) * np.log(x)


MODELS = {
    "exponential": (_exp_logpdf, _exp_logsf),
    "weibull": (_weib_logpdf, _weib_logsf),
    "lognormal": (_ln_logpdf, _ln_logsf),
    "pareto": (_par_logpdf, _par_logsf),
}


def _start(name, x, u):
    lx = np.log(x)
    if name == "exponential":
        return [np.log(1.0 / max(np.mean(x - u), 1e-9))]
    if name == "weibull":
        return [np.log(np.mean(x)), np.log(0.7)]
    if name == "lognormal":
        return [np.mean(lx), np.log(max(np.std(lx), 0.05))]
    if name == "pareto":
        return [np.log(1.0 / max(np.mean(lx - np.log(u)), 1e-9))]


def truncated_loglik_pointwise(name, th, x, u):
    lp, ls = MODELS[name]
    return lp(x, th) - ls(u, th)


def fit_truncated(name, x, u, n_cens=0, c_cens=None):
    """MLE of the model on exceedances x > u. Returns dict(theta, params, loglik, pointwise).

    n_cens observations right-censored at c_cens (> u) enter via log P(S > c_cens | S > u).
    """
    x = np.asarray(x, dtype=np.float64)
    lsf = MODELS[name][1]

    def f(th):
        v = -np.sum(truncated_loglik_pointwise(name, th, x, u))
        if n_cens:
            v -= n_cens * (lsf(c_cens, th) - lsf(u, th))
        return v
    best = None
    starts = [_start(name, x, u)]
    if name == "weibull":
        starts += [[np.log(u), np.log(b)] for b in (0.3, 0.5, 1.0, 2.0)]
    if name == "lognormal":
        starts += [[np.log(u) - k, np.log(s)] for k in (0, 2, 5) for s in (0.3, 1.0, 2.0)]
    for s0 in starts:
        with np.errstate(all="ignore"):
            r = optimize.minimize(f, s0, method="Nelder-Mead",
                                  options=dict(maxiter=20000, xatol=1e-9, fatol=1e-9))
        if np.isfinite(r.fun) and (best is None or r.fun < best.fun):
            best = r
    th = best.x
    pw = truncated_loglik_pointwise(name, th, x, u)
    return dict(theta=th, params=named_params(name, th), loglik=float(pw.sum()), pointwise=pw, k=len(th))


def named_params(name, th):
    if name == "exponential":
        return dict(theta=float(np.exp(th[0])))
    if name == "weibull":
        return dict(lam=float(np.exp(th[0])), beta=float(np.exp(th[1])))
    if name == "lognormal":
        return dict(mu=float(th[0]), sigma=float(np.exp(th[1])))
    if name == "pareto":
        return dict(alpha=float(np.exp(th[0])))


def vuong(fa, fb):
    """Vuong normalized LR test between two fitted models (non-nested form).

    Returns (R, z, p_two_sided). R > 0 favours model a.
    """
    d = fa["pointwise"] - fb["pointwise"]
    n = len(d)
    R = d.sum()
    sd = d.std(ddof=1)
    if sd == 0:
        return R, 0.0, 1.0
    z = R / (np.sqrt(n) * sd)
    return float(R), float(z), float(2 * stats.norm.sf(abs(z)))


def fit_gpd(x, u):
    """Generalized Pareto fit to excesses x-u. Returns (xi, scale)."""
    y = np.asarray(x, dtype=np.float64) - u
    xi, loc, scale = stats.genpareto.fit(y, floc=0)
    return float(xi), float(scale)


def gpd_xi_ci(x, u, nboot=200, rng=None):
    rng = np.random.default_rng(0) if rng is None else rng
    y = np.asarray(x, dtype=np.float64)
    xs = []
    for _ in range(nboot):
        b = rng.choice(y, size=len(y), replace=True)
        xs.append(fit_gpd(b, u)[0])
    return np.percentile(xs, [2.5, 97.5])


# ------------------------------------------------------------------ tail model evaluation
def model_sf(name, params, x, u, pu):
    """Unconditional survival P(S > x) = P(S > u) * P(S > x | S > u) under a fitted tail model."""
    x = np.asarray(x, dtype=np.float64)
    if name == "exponential":
        r = np.exp(-params["theta"] * (x - u))
    elif name == "weibull":
        lam, b = params["lam"], params["beta"]
        r = np.exp(-((x / lam) ** b) + (u / lam) ** b)
    elif name == "lognormal":
        mu, s = params["mu"], params["sigma"]
        r = np.exp(special.log_ndtr(-(np.log(x) - mu) / s) - special.log_ndtr(-(np.log(u) - mu) / s))
    elif name == "pareto":
        r = (x / u) ** (-params["alpha"])
    elif name == "gpd":
        r = stats.genpareto.sf(x - u, params["xi"], scale=params["scale"])
    return pu * np.minimum(r, 1.0)


def model_isf(name, params, eps, u, pu):
    """Inverse of model_sf: x with P(S > x) = eps (eps < pu)."""
    q = np.asarray(eps, dtype=np.float64) / pu  # conditional exceedance prob
    lq = np.log(q)
    if name == "exponential":
        return u - lq / params["theta"]
    if name == "weibull":
        lam, b = params["lam"], params["beta"]
        return lam * ((u / lam) ** b - lq) ** (1 / b)
    if name == "lognormal":
        mu, s = params["mu"], params["sigma"]
        base = special.log_ndtr(-(np.log(u) - mu) / s)
        z = -special.ndtri(np.exp(lq + base))
        return np.exp(mu + s * z)
    if name == "pareto":
        return u * q ** (-1 / params["alpha"])
    if name == "gpd":
        return u + stats.genpareto.isf(q, params["xi"], scale=params["scale"])


def analyze_tail(S, u=None, tail_frac=1e-2, min_exceed=200, censored_at=None):
    """Run the full tail-class pipeline on samples S. Returns a dict of results.

    u defaults to the empirical (1 - tail_frac) quantile (raised if needed so there are at
    least min_exceed exceedances). censored_at: values >= this are right-censored (Relay leg limit);
    they are excluded from fits (reported separately) — callers should ensure the censored
    fraction is far below the fitted tail probabilities.
    """
    S = np.asarray(S, dtype=np.float64)
    n = len(S)
    if u is None:
        u = quantile_tail(S, tail_frac)
    x = S[S > u] + 0.5  # continuity correction
    ncens = 0
    if censored_at is not None:
        ncens = int((S >= censored_at).sum())
        x = x[x < censored_at]
    pu = (S > u).mean()
    out = dict(n=n, u=float(u), pu=float(pu), n_exceed=int(len(x)), n_censored=ncens, fits={}, vuong={})
    if len(x) < 20:
        return out
    for name in MODELS:
        out["fits"][name] = fit_truncated(name, x, u + 0.5, n_cens=ncens, c_cens=censored_at)
    try:
        xi, sc = fit_gpd(x, u + 0.5)
        ci = gpd_xi_ci(x, u + 0.5, nboot=100)
        out["gpd"] = dict(xi=xi, scale=sc, xi_ci=[float(ci[0]), float(ci[1])])
    except Exception as e:  # pragma: no cover
        out["gpd"] = dict(error=str(e))
    names = list(MODELS)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            out["vuong"][f"{a}|{b}"] = vuong(out["fits"][a], out["fits"][b])
    # AIC ranking
    aic = {k: 2 * f["k"] - 2 * f["loglik"] for k, f in out["fits"].items()}
    out["aic"] = aic
    out["best_aic"] = min(aic, key=aic.get)
    for f in out["fits"].values():
        f.pop("pointwise", None)
        f["theta"] = [float(t) for t in f["theta"]]
    return out
