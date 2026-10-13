"""Tail characterization of decoder work distributions (Experiment 2).

Peaks over threshold: the excesses Y = S - u of the samples above a high threshold u are fitted by
maximum likelihood to a generalized Pareto distribution GPD(xi, sigma). The shape xi is the tail
parameter: xi > 0 is a heavy (power-law-type) tail with P(S > x) ~ x^(-1/xi); xi = 0 is the
exponential tail; xi < 0 is a bounded tail. Only two likelihood-ratio procedures are used:

  * LR test of H0: xi = 0 (exponential excesses) against free xi, chi^2_1 under H0;
  * profile-likelihood 95% confidence interval for xi (the inverted LR test).

Work counts are integers; samples above u are shifted by +0.5 (continuity correction), which is
immaterial when u is in the tens or more. Right-censored samples (Relay-BP leg limit) enter the
likelihood through the GPD survival function.
"""
import numpy as np
from scipy import optimize, stats

CHI2_95 = stats.chi2.ppf(0.95, 1)
_XI_EPS = 1e-9


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


# ------------------------------------------------------------------ GPD likelihood
def gpd_loglik(xi, sigma, y, n_cens=0, y_cens=None):
    """Log-likelihood of excesses y (and n_cens excesses right-censored at y_cens) under GPD(xi, sigma)."""
    if sigma <= 0:
        return -np.inf
    if abs(xi) < _XI_EPS:
        ll = -len(y) * np.log(sigma) - y.sum() / sigma
        if n_cens:
            ll += n_cens * (-y_cens / sigma)
        return ll
    t = 1 + xi * y / sigma
    if np.any(t <= 0) or (n_cens and 1 + xi * y_cens / sigma <= 0):
        return -np.inf
    ll = -len(y) * np.log(sigma) - (1 + 1 / xi) * np.log1p(xi * y / sigma).sum()
    if n_cens:
        ll += n_cens * (-np.log1p(xi * y_cens / sigma) / xi)
    return ll


def _profile(xi, y, n_cens, y_cens):
    """max over sigma of the log-likelihood at fixed xi; returns (loglik, sigma)."""
    if abs(xi) < _XI_EPS:
        s = (y.sum() + n_cens * (y_cens or 0.0)) / len(y)  # closed-form exponential MLE
        return gpd_loglik(0.0, s, y, n_cens, y_cens), s
    ymax = max(y.max(), y_cens if n_cens else 0.0)
    lo = np.log(max(-xi * ymax, 1e-12) * (1 + 1e-9)) if xi < 0 else np.log(y.mean()) - 10
    hi = np.log(y.mean()) + 10
    r = optimize.minimize_scalar(lambda ls: -gpd_loglik(xi, np.exp(ls), y, n_cens, y_cens),
                                 bounds=(lo, max(hi, lo + 1)), method="bounded",
                                 options=dict(xatol=1e-10))
    return -r.fun, float(np.exp(r.x))


def fit_gpd(x, u, n_cens=0, c_cens=None):
    """MLE of GPD(xi, sigma) on samples x > u, plus the LR test of xi = 0 and the profile-likelihood CI.

    Returns dict(xi, scale, loglik, xi_ci=[lo, hi], lr_exp=dict(stat, p)).
    """
    y = np.asarray(x, dtype=np.float64) - u
    y_cens = (c_cens - u) if n_cens else None
    # starting points: unconstrained scipy fit and a grid of shapes
    starts = [0.0, 0.1, 0.3, 0.6, 1.0, -0.1]
    try:
        starts.append(float(stats.genpareto.fit(y, floc=0)[0]))
    except Exception:  # pragma: no cover
        pass
    prof = lambda xi: _profile(xi, y, n_cens, y_cens)[0]
    best_xi = max(starts, key=prof)
    r = optimize.minimize_scalar(lambda xi: -prof(xi), bracket=(best_xi - 0.05, best_xi + 0.05),
                                 options=dict(xtol=1e-8))
    xi_hat = float(r.x) if -r.fun >= prof(best_xi) else best_xi
    ll_hat, sigma_hat = _profile(xi_hat, y, n_cens, y_cens)
    ll_exp = prof(0.0)
    stat = max(2 * (ll_hat - ll_exp), 0.0)

    def ci_end(direction):
        step, xi0 = 0.02, xi_hat
        g = lambda xi: ll_hat - prof(xi) - CHI2_95 / 2
        for _ in range(200):
            xi1 = xi0 + direction * step
            if g(xi1) > 0:
                return float(optimize.brentq(g, min(xi0, xi1), max(xi0, xi1), xtol=1e-6))
            xi0, step = xi1, step * 1.5
        return float(direction * np.inf)

    return dict(xi=xi_hat, scale=sigma_hat, loglik=float(ll_hat), xi_ci=[ci_end(-1), ci_end(+1)],
                lr_exp=dict(stat=float(stat), p=float(stats.chi2.sf(stat, 1))))


# ------------------------------------------------------------------ fitted-tail evaluation
def gpd_sf(params, x, u, pu):
    """Unconditional P(S > x) = P(S > u) * P(S > x | S > u) under the fitted GPD tail."""
    x = np.asarray(x, dtype=np.float64)
    return pu * stats.genpareto.sf(x - u, params["xi"], scale=params["scale"])


def gpd_isf(params, eps, u, pu):
    """Inverse of gpd_sf: x with P(S > x) = eps (eps < pu)."""
    q = np.asarray(eps, dtype=np.float64) / pu
    return u + stats.genpareto.isf(q, params["xi"], scale=params["scale"])


def analyze_tail(S, u=None, tail_frac=1e-2, censored_at=None):
    """GPD tail fit of samples S above the (1 - tail_frac) empirical quantile (or the given u).

    censored_at: samples >= this value are right-censored (Relay-BP leg limit) and enter the
    likelihood through the GPD survival function.
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
    out = dict(n=n, u=float(u), pu=float(pu), n_exceed=int(len(x)), n_censored=ncens)
    if len(x) < 20:
        return out
    out["gpd"] = fit_gpd(x, u + 0.5, n_cens=ncens, c_cens=censored_at)
    return out
