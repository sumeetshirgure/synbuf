"""Run all Experiment 1 / Experiment 2 analyses on collected data.

Outputs: results/*.json, results/*.csv and figures/*.png.
Usage: PYTHONPATH=. python3 scripts/analyze.py [--skip-lindley]
"""
import argparse
import warnings
warnings.filterwarnings("ignore")
import csv
import glob
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from synbuf import buffer as buf  # noqa: E402
from synbuf import exp1, tails  # noqa: E402
from synbuf.circuits import BB_CODES  # noqa: E402
from synbuf.run_relay import load as load_relay  # noqa: E402

RES, FIG = "results", "figures"
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
EPS_GRID = np.logspace(-6, -15, 10)
MODEL_NAMES = ["exponential", "lognormal", "pareto"]
COLORS = {"exponential": "#1f77b4", "lognormal": "#d62728", "pareto": "#9467bd",
          "gpd": "#ff7f0e"}
plt.rcParams.update({"font.size": 9, "figure.dpi": 130, "axes.grid": True, "grid.alpha": 0.3})


class NpEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def dump(obj, name):
    with open(os.path.join(RES, name), "w") as f:
        json.dump(obj, f, indent=1, cls=NpEncoder)


# ------------------------------------------------------------------ data loading
def load_mwpm():
    pts = {}
    for f in sorted(glob.glob("data/mwpm/surface_d*_p*.npz")):
        m = re.search(r"d(\d+)_p([0-9.e-]+)\.npz", f)
        d, p = int(m.group(1)), float(m.group(2))
        z = np.load(f)
        pts[(d, p)] = dict(S=z["S"], pq=z["pq"], fail=z["fail"],
                           label=f"surface d={d}", d=d, p=p, decoder="MWPM")
    return pts


def load_bb(tag=""):
    pts = {}
    for code in BB_CODES:
        for dd in sorted(glob.glob(f"data/relay/bb{code}_p*{tag}")):
            m = re.search(r"_p([0-9.e-]+)" + re.escape(tag) + "$", dd)
            if not m:
                continue
            p = float(m.group(1))
            r = load_relay(dd)
            if r is None:
                continue
            meta = r["meta"]
            cens_at = meta["pre_iter"] + meta["num_sets"] * meta["set_max_iter"]
            pts[(code, p)] = dict(S=r["it_first"], S_total=r["it_total"], fail=r["fail"],
                                  fail_first=r["fail_first"], conv=r["conv"], legs=r["legs"],
                                  meta=meta, cens_at=cens_at,
                                  label=f"[[{code.replace('_', ',')}]]", d=BB_CODES[code].d, p=p,
                                  decoder="Relay-BP", code=code)
    return pts


# ------------------------------------------------------------------ Experiment 1
def _order(k):
    """Sort points by distance, then block length (BB codes) / key."""
    if isinstance(k[0], str):
        c = BB_CODES[k[0]]
        return (c.d, c.n)
    return (k[0], 0)


def experiment1(pts, ref_key, family, extra_caps=None, xkey="d"):
    S_ref = pts[ref_key]["S"]
    caps = {f"c{j}": float(tails.quantile_tail(S_ref, 10.0 ** -j)) for j in (3, 4, 5)}
    not_calibratable = {}
    cens = pts[ref_key].get("cens_at")
    n_ref = len(S_ref)
    for j in (3, 4, 5):
        name = f"c{j}"
        # a cap needs >= 10 shots above it at the reference point and must lie below any censoring limit
        if (cens is not None and caps[name] >= cens) or 10.0 ** -j * n_ref < 10:
            not_calibratable[name] = dict(
                reason="P(S > c) = 1e-%d is beyond the leg limit or the sample size at the reference point" % j,
                censored_fraction=float((S_ref >= cens).mean()) if cens else 0.0, n_ref=n_ref)
            caps.pop(name)
    if extra_caps:
        caps.update(extra_caps)
    rows = []
    for key, P in sorted(pts.items(), key=lambda kv: (kv[1]["p"], kv[1]["d"], str(kv[0]))):
        for cname, c in caps.items():
            st = exp1.cap_stats(P["S"], P["fail"], c)
            rows.append(dict(family=family, point=P["label"], key=str(key), d=P["d"], p=P["p"], cap=cname,
                             **{k: v for k, v in st.items()}))
    dump(dict(reference=str(ref_key), caps=caps, not_calibratable=not_calibratable, rows=rows), f"exp1_{family}.json")
    with open(os.path.join(RES, f"exp1_{family}.csv"), "w", newline="") as f:
        cols = ["point", "d", "p", "cap", "c", "n", "pL0", "p_slow", "PL", "dPL", "p_fail_given_slow",
                "k_slow", "k_fail_slow", "k_fail"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return caps, rows


def fig_ler_vs_d(pts, caps, rows, family, title, ps=None, xlabels=None):
    ps = sorted({P["p"] for P in pts.values()}) if ps is None else ps
    fig, axes = plt.subplots(1, len(ps), figsize=(3.3 * len(ps), 3.1), sharey=True, squeeze=False)
    capcol = {"c3": "#d62728", "c4": "#ff7f0e", "c5": "#2ca02c", "deployed": "#9467bd"}
    for ax, p in zip(axes[0], ps):
        keys = sorted([k for k, P in pts.items() if P["p"] == p], key=_order)
        if not keys:
            continue
        xs = np.arange(len(keys)) if xlabels else np.array([pts[k]["d"] for k in keys])
        pL0 = np.array([pts[k]["fail"].mean() for k in keys])
        lo = np.array([exp1.wilson(int(pts[k]["fail"].sum()), len(pts[k]["fail"]))[0] for k in keys])
        hi = np.array([exp1.wilson(int(pts[k]["fail"].sum()), len(pts[k]["fail"]))[1] for k in keys])
        ax.errorbar(xs, np.where(pL0 > 0, pL0, np.nan), yerr=[np.maximum(pL0 - lo, 0), np.maximum(hi - pL0, 0)], fmt="k--o", ms=3,
                    label=r"uncapped $p_L^{(0)}$", capsize=2)
        for cname, c in caps.items():
            r = [next(x for x in rows if x["key"] == str(k) and x["cap"] == cname) for k in keys]
            y = np.array([x["PL"] for x in r])
            yl = np.array([x["PL_ci"][0] for x in r])
            yh = np.array([x["PL_ci"][1] for x in r])
            ax.errorbar(xs, y, yerr=[np.maximum(y - yl, 0), np.maximum(yh - y, 0)], fmt="-o", ms=3, capsize=2, color=capcol.get(cname),
                        label=f"$P_L(c)$, {cname}={c:.0f}")
        ax.set_yscale("log")
        ax.set_title(f"{title}, p={p:g}")
        if xlabels:
            ax.set_xticks(xs)
            ax.set_xticklabels([xlabels(k) for k in keys], rotation=30, fontsize=7)
        else:
            ax.set_xlabel("code distance d")
    axes[0][0].set_ylabel("logical error rate per window")
    axes[0][-1].legend(fontsize=6.5, loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"exp1_ler_vs_d_{family}.png"))
    plt.close(fig)


def fig_PL_vs_c(pts, keys, family, title):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for k in keys:
        P = pts[k]
        S = P["S"]
        cs = np.unique(np.quantile(S, np.linspace(0.5, 1, 400)))
        cs = np.unique(np.r_[cs, np.sort(S)[-2000:]])
        cc = exp1.cap_curve(S, P["fail"], cs)
        l, = axes[0].loglog(cc["c"], cc["PL"], label=P["label"] + f" p={P['p']:g}")
        axes[0].axhline(cc["pL0"], color=l.get_color(), ls=":", lw=1)
        ok = cc["k_slow"] >= 5
        axes[1].semilogx(cc["c"][ok], cc["p_fail_given_slow"][ok], color=l.get_color())
    axes[0].set_xlabel("cap c (work units)")
    axes[0].set_ylabel(r"$P_L(c)$ (dotted: $p_L^{(0)}$)")
    axes[0].set_title(title + r": $P_L(c)$ vs $c$")
    axes[1].set_xlabel("cap c (work units)")
    axes[1].set_ylabel(r"$P(F_\infty \mid S>c)$")
    axes[1].set_title("are slow shots the failing shots?")
    axes[1].set_ylim(0, 1)
    axes[0].legend(fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"exp1_PL_vs_c_{family}.png"))
    plt.close(fig)


def scaled_cap_rows(pts, family, ks=(3, 5, 8, 12)):
    """Caps that scale with the typical work: c = k * median(S) at each point."""
    rows = []
    for key, P in sorted(pts.items(), key=lambda kv: (kv[1]["p"], kv[1]["d"], str(kv[0]))):
        med = float(np.median(P["S"]))
        for k in ks:
            st = exp1.cap_stats(P["S"], P["fail"], k * med)
            rows.append(dict(family=family, point=P["label"], key=str(key), d=P["d"], p=P["p"], k=k,
                             median=med, **st))
    dump(rows, f"exp1_scaled_caps_{family}.json")
    return rows


def fig_scaled_caps(pts, rows, family, title, xlabels=None):
    ps = sorted({P["p"] for P in pts.values()})
    fig, axes = plt.subplots(1, len(ps), figsize=(3.3 * len(ps), 3.1), sharey=True, squeeze=False)
    for ax, p in zip(axes[0], ps):
        keys = sorted([k for k, P in pts.items() if P["p"] == p], key=_order)
        xs = np.arange(len(keys)) if xlabels else np.array([pts[k]["d"] for k in keys])
        pL0 = np.array([pts[k]["fail"].mean() for k in keys])
        ax.plot(xs, np.where(pL0 > 0, pL0, np.nan), "k--o", ms=3, label=r"uncapped $p_L^{(0)}$")
        for kk in sorted({r["k"] for r in rows}):
            y = [next(r for r in rows if r["key"] == str(k) and r["k"] == kk)["PL"] for k in keys]
            ax.plot(xs, y, "-o", ms=3, label=f"$P_L(c)$, c={kk}x median(S)")
        ax.set_yscale("log")
        ax.set_title(f"{title}, p={p:g}")
        if xlabels:
            ax.set_xticks(xs)
            ax.set_xticklabels([xlabels(k) for k in keys], rotation=30, fontsize=7)
        else:
            ax.set_xlabel("code distance d")
    axes[0][0].set_ylabel("logical error rate per window")
    axes[0][-1].legend(fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"exp1_scaled_caps_{family}.png"))
    plt.close(fig)


# ------------------------------------------------------------------ Experiment 2
def tail_analysis(P, tail_frac):
    cens = P.get("cens_at")
    out = tails.analyze_tail(P["S"], tail_frac=tail_frac, censored_at=cens)
    S = P["S"]
    out["mean"] = float(S.mean())
    out["median"] = float(np.median(S))
    out["std"] = float(S.std())
    out["max"] = float(S.max())
    out["q"] = {f"1e-{j}": float(tails.quantile_tail(S, 10.0 ** -j)) for j in range(1, 8)
                if 10.0 ** -j * len(S) >= 1}
    out["censored_total"] = int((S >= cens).sum()) if cens else 0
    return out


def _empirical_tail_integral(S, u):
    """Return g(x) = int_x^u Fbar_emp(y) dy for x <= u (exact for the empirical step function)."""
    x, sf = tails.survival(S)
    keep = x < u
    x, sf = x[keep], sf[keep]
    grid = np.r_[x, u]
    seg = sf * np.diff(grid)  # Fbar is constant (= sf_i) on [x_i, x_{i+1})
    cum = np.r_[np.cumsum(seg[::-1])[::-1], 0.0]  # integral from grid_i to u

    def g(t):
        if t >= u:
            return 0.0
        i = np.searchsorted(grid, t, side="right") - 1
        if i < 0:
            return cum[0] + (grid[0] - t)  # Fbar = 1 below the minimum
        return cum[i + 1] + sf[i] * (grid[i + 1] - t)
    return g


def _buffer_rows(sf_tail, S, u, rho, ms, isf=None):
    S = np.asarray(S, dtype=np.float64)
    ES, med = S.mean(), np.median(S)
    g = _empirical_tail_integral(S, u)
    f = lambda t: sf_tail(np.exp(t)) * np.exp(t)
    from scipy import integrate
    tail_int = 0.0
    lo = np.log(u)
    for k in range(40):  # integrate in decades until negligible
        val, _ = integrate.quad(f, lo + k, lo + k + 1, limit=200, epsabs=0, epsrel=1e-7)
        tail_int += val
        if val < 1e-30 * max(tail_int, 1e-300) or (k > 3 and val < 1e-6 * tail_int):
            break
    if not np.isfinite(tail_int):
        return None

    def FI(x):
        if x < u:
            return (g(x) + tail_int) / ES
        val = 0.0
        lx = np.log(x)
        for k in range(40):
            v, _ = integrate.quad(f, lx + k, lx + k + 1, limit=200, epsabs=0, epsrel=1e-7)
            val += v
            if k > 3 and v < 1e-6 * val:
                break
        return val / ES

    x_emp, sf_emp = tails.survival(S)
    rows = []
    for eps in EPS_GRID:
        r = dict(eps=float(eps))
        if isf is not None:
            r["x_quantile"] = float(isf(eps))
        for m in ms:
            tgt = eps * (1 - rho) / rho if m == 1 else eps ** (1.0 / m)
            r[f"x_m{m}"] = buf.solve_x(FI, tgt, x_lo=max(S.min(), 1.0), x_hi=1e40)
        for key in [k for k in r if k.startswith("x_")]:
            r[key.replace("x_", "med_")] = r[key] / med
            r[key.replace("x_", "win_")] = rho * r[key] / ES
        rows.append(r)
    return rows


def buffer_table(P, ta, rho=0.5, ms=(1, 2, 3)):
    """Buffer requirements from each fitted tail model, spliced onto the empirical bulk below u.

    x_m1: D/G/1 (sliding window), solves rho/(1-rho) * Fbar_I(x) = eps   (Veraverbeke)
    x_mK: D/G/k with K spare servers, solves Fbar_I(x)^K = eps           (Foss-Korshunov, leading order)
    x_quantile: plain Fbar^{-1}(eps) (EXPERIMENTS.md table leading order)
    """
    S = P["S"].astype(np.float64)
    u, pu = ta["u"] + 0.5, ta["pu"]
    out = {}
    for name in list(ta["fits"]) + (["gpd"] if "xi" in ta.get("gpd", {}) else []):
        params = ta["fits"][name]["params"] if name != "gpd" else ta["gpd"]
        if name == "pareto" and params["alpha"] <= 1:
            out[name] = dict(note="alpha<=1: infinite mean, integrated tail diverges")
            continue
        if name == "gpd" and params["xi"] >= 1:
            out[name] = dict(note="xi>=1: infinite mean")
            continue
        sf = lambda x, name=name, params=params: float(tails.model_sf(name, params, x, u, pu))
        isf = lambda e, name=name, params=params: tails.model_isf(name, params, e, u, pu)
        rows = _buffer_rows(sf, S, u, rho, ms, isf=isf)
        out[name] = rows if rows is not None else dict(note="integrated tail diverges")
    return out


def fig_diagnostics(P, ta, fname, title):
    S = P["S"].astype(np.float64)
    x, sf = tails.survival(S)
    ok = sf > 0
    x, sf = x[ok], sf[ok]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.0))
    ax = axes[0]
    ax.semilogy(x, sf, "k.", ms=2)
    ax.set_xlabel("S")
    ax.set_ylabel(r"$\bar F(S)$")
    ax.set_title(r"$\log\bar F$ vs $x$ (exponential: line)")
    ax = axes[1]
    ax.loglog(x, sf, "k.", ms=2, label="empirical")
    u, pu = ta["u"] + 0.5, ta["pu"]
    xx = np.logspace(np.log10(u), np.log10(x.max() * 30), 200)
    for name, f in ta["fits"].items():
        ax.loglog(xx, tails.model_sf(name, f["params"], xx, u, pu), color=COLORS[name], lw=1, label=name)
    ax.set_ylim(max(sf.min() / 30, 1e-12), 1.2)
    ax.set_xlabel("S")
    ax.set_title(r"$\log\bar F$ vs $\log x$ (Pareto: line) + tail fits")
    ax.legend(fontsize=6)
    ax = axes[2]
    lS = np.log(np.maximum(S, 1))
    sub = np.sort(lS)[:: max(1, len(lS) // 20000)]
    tailpart = np.sort(lS)[-2000:]
    qq = np.r_[sub, tailpart]
    n = len(lS)
    ranks = np.r_[np.arange(0, n, max(1, n // 20000))[: len(sub)], np.arange(n - 2000, n)]
    theo = stats.norm.ppf((ranks + 0.5) / n)
    ax.plot(theo, qq, "k.", ms=2)
    ax.set_xlabel("normal quantile")
    ax.set_ylabel(r"$\log S$")
    ax.set_title(r"normal Q-Q of $\log S$ (lognormal: line)")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname))
    plt.close(fig)


def fig_extrapolation(P, ta, fname, title):
    S = P["S"].astype(np.float64)
    x, sf = tails.survival(S)
    ok = sf > 0
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    ax.loglog(x[ok], sf[ok], "k.", ms=2, label="empirical")
    u, pu = ta["u"] + 0.5, ta["pu"]
    xmax = max(x.max() * 30, 10 * u)
    for name, f in ta["fits"].items():
        try:
            xe = tails.model_isf(name, f["params"], 1e-15, u, pu)
            xmax = max(xmax, min(xe, 1e9))
        except Exception:
            pass
    xx = np.logspace(np.log10(u), np.log10(xmax), 300)
    for name, f in ta["fits"].items():
        ax.loglog(xx, tails.model_sf(name, f["params"], xx, u, pu), color=COLORS[name], lw=1, label=name + " fit")
    ax.set_ylim(1e-16, 1.5)
    ax.axhline(1e-12, color="gray", lw=0.6, ls=":")
    ax.set_xlabel("S (work units)")
    ax.set_ylabel(r"$P(S>x)$")
    ax.set_title(title, fontsize=8)
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-lindley", action="store_true")
    a = ap.parse_args()
    summary = {}

    # ---------------- MWPM
    mw = load_mwpm()
    bb = load_bb()
    print("MWPM points:", len(mw), " BB points:", len(bb))

    if mw:
        ref = (9, 0.003) if (9, 0.003) in mw else sorted(mw)[len(mw) // 2]
        caps, rows = experiment1(mw, ref, "mwpm")
        summary["mwpm_caps"] = dict(reference=ref, caps=caps)
        fig_ler_vs_d(mw, caps, rows, "mwpm", "MWPM surface code")
        keys = [k for k in [(9, 0.003), (11, 0.003), (13, 0.003), (15, 0.003)] if k in mw]
        fig_PL_vs_c(mw, keys, "mwpm", "MWPM")
        fig_scaled_caps(mw, scaled_cap_rows(mw, "mwpm"), "mwpm", "MWPM, scaled cap")

    if bb:
        ref = ("144_12_12", 0.003) if ("144_12_12", 0.003) in bb else sorted(bb)[0]
        dep = bb[ref]["meta"]["pre_iter"] + 300 * bb[ref]["meta"]["set_max_iter"]
        caps, rows = experiment1(bb, ref, "relay", extra_caps={"deployed": float(dep)})
        summary["relay_caps"] = dict(reference=ref, caps=caps)
        fig_ler_vs_d(bb, caps, rows, "relay", "Relay-BP BB",
                     xlabels=lambda k: f"{BB_CODES[k[0]].n},{BB_CODES[k[0]].k},{BB_CODES[k[0]].d}")
        keys = [k for k in sorted(bb, key=lambda k: BB_CODES[k[0]].n) if k[1] == 0.003]
        fig_PL_vs_c(bb, keys, "relay", "Relay-BP")
        fig_scaled_caps(bb, scaled_cap_rows(bb, "relay"), "relay", "Relay-BP, scaled cap",
                        xlabels=lambda k: f"{BB_CODES[k[0]].n},{BB_CODES[k[0]].k},{BB_CODES[k[0]].d}")

    # ---------------- Experiment 2: tails
    tail_rows = []
    allta = {}
    for fam, pts, frac in [("mwpm", mw, 1e-3), ("relay", bb, 1e-2)]:
        for key, P in sorted(pts.items(), key=lambda kv: (kv[1]["p"], kv[1]["d"], str(kv[0]))):
            if fam == "mwpm" and P["d"] < 5:
                continue
            fr = frac if len(P["S"]) * frac >= 2000 else 2000 / len(P["S"])
            ta = tail_analysis(P, fr)
            ta2 = tails.analyze_tail(P["S"], tail_frac=fr / 10, censored_at=P.get("cens_at"))
            ta["robust_threshold"] = dict(u=ta2["u"], n_exceed=ta2["n_exceed"], best_aic=ta2.get("best_aic"),
                                          gpd=ta2.get("gpd"),
                                          params={k: v["params"] for k, v in ta2.get("fits", {}).items()})
            allta[(fam, key)] = ta
            row = dict(family=fam, point=P["label"], d=P["d"], p=P["p"], n=ta["n"], mean=ta["mean"],
                       median=ta["median"], u=ta["u"], n_exceed=ta["n_exceed"], best_aic=ta.get("best_aic"),
                       xi=ta.get("gpd", {}).get("xi"), xi_lo=ta.get("gpd", {}).get("xi_ci", [None, None])[0],
                       xi_hi=ta.get("gpd", {}).get("xi_ci", [None, None])[1],
                       xi_high_u=(ta2.get("gpd") or {}).get("xi"), best_aic_high_u=ta2.get("best_aic"),
                       censored=ta["censored_total"])
            for name in MODEL_NAMES:
                if name in ta["fits"]:
                    for pk, pv in ta["fits"][name]["params"].items():
                        row[f"{name}_{pk}"] = pv
                    row[f"{name}_aic"] = ta["aic"][name]
            for pair, (R, z, pv) in ta["vuong"].items():
                row[f"vuong_{pair}_z"] = z
                row[f"vuong_{pair}_p"] = pv
            tail_rows.append(row)
            print(f"[tail] {fam} {P['label']} p={P['p']:g}: best={row['best_aic']} xi={row['xi']:.3f} "
                  f"(high-u xi={row['xi_high_u']})", flush=True)
    dump(tail_rows, "exp2_tail_fits.json")
    if tail_rows:
        cols = sorted({k for r in tail_rows for k in r}, key=lambda c: (c not in ("family", "point", "d", "p"), c))
        with open(os.path.join(RES, "exp2_tail_fits.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(tail_rows)

    # ---------------- diagnostics / extrapolation / buffers for selected points
    sel = []
    for k in [(15, 0.003), (11, 0.003), (15, 0.001), (15, 0.0005)]:
        if k in mw:
            sel.append(("mwpm", k, mw[k], f"MWPM surface d={k[0]} p={k[1]:g}"))
    for k in [("144_12_12", 0.003), ("108_8_10", 0.003), ("90_8_10", 0.003), ("72_12_6", 0.003), ("144_12_12", 0.001)]:
        if k in bb:
            sel.append(("relay", k, bb[k], f"Relay-BP [[{k[0].replace('_', ',')}]] p={k[1]:g}"))
    buffers = {}
    for fam, key, P, title in sel:
        ta = allta[(fam, key)]
        tag = f"{fam}_{'_'.join(str(x) for x in key)}"
        fig_diagnostics(P, ta, f"exp2_diag_{tag}.png", title)
        fig_extrapolation(P, ta, f"exp2_extrap_{tag}.png", title)
        buffers[tag] = dict(title=title, mean=float(P["S"].mean()), median=float(np.median(P["S"])),
                            best=ta.get("best_aic"), table=buffer_table(P, ta))
    dump(buffers, "exp2_buffers.json")

    nconv5_analysis(bb)
    fig_param_trends(tail_rows)
    fig_buffers(buffers)

    # ---------------- Lindley check (D/G/1) for selected points
    if not a.skip_lindley:
        lind = {}
        for fam, key, P, title in sel[:1] + [s for s in sel if s[0] == "relay"][:1]:
            S = P["S"].astype(np.float64)
            thr_med = np.array([5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
            thr = thr_med * np.median(S)
            res = {}
            for rho in (0.5, 0.8):
                pr, aa = buf.lindley_overflow(S, rho, thr, n_windows=200_000_000)
                ta = allta[(fam, key)]
                u, pu = ta["u"] + 0.5, ta["pu"]
                best = ta["best_aic"]
                sf = lambda x: tails.model_sf(best, ta["fits"][best]["params"], x, u, pu)
                emp_sf = lambda x: np.mean(S > x)
                FI_emp = [np.mean(np.maximum(S - x, 0)) / S.mean() for x in thr]
                res[str(rho)] = dict(a=aa, thresholds_med=thr_med.tolist(), thresholds=thr.tolist(), P_overflow=pr.tolist(),
                                     veraverbeke_empirical=[rho / (1 - rho) * f for f in FI_emp])
            lind[f"{fam}_{'_'.join(str(x) for x in key)}"] = res
            print("[lindley]", key, res, flush=True)
        dump(lind, "exp2_lindley.json")
        fig_lindley(lind)
    dump({k: (v if not isinstance(v, dict) else {kk: (str(vv) if isinstance(vv, tuple) else vv) for kk, vv in v.items()})
          for k, v in summary.items()}, "summary.json")


def nconv5_analysis(bb1):
    """Time to first solution (S_1) vs time to collect 5 solutions (S_5), same shots, stop_nconv=5."""
    b5 = load_bb("_nconv5")
    if not b5:
        return
    out = {}
    fig, axes = plt.subplots(1, len(b5), figsize=(4.0 * len(b5), 3.2), squeeze=False)
    for ax, (key, P) in zip(axes[0], sorted(b5.items())):
        S1, S5 = P["S"].astype(float), P["S_total"].astype(float)
        r = dict(n=len(S1), ler_best_of_5=float(P["fail"].mean()), ler_first=float(P["fail_first"].mean()),
                 censored=int((~P["conv"]).sum()), mean_S1=float(S1.mean()), mean_S5=float(S5.mean()),
                 median_S1=float(np.median(S1)), median_S5=float(np.median(S5)),
                 q_S1={f"1e-{j}": float(tails.quantile_tail(S1, 10.0 ** -j)) for j in (1, 2, 3, 4) if 10.0 ** -j * len(S1) >= 5},
                 q_S5={f"1e-{j}": float(tails.quantile_tail(S5, 10.0 ** -j)) for j in (1, 2, 3, 4) if 10.0 ** -j * len(S5) >= 5})
        if key in bb1:
            r["ler_nconv1_main_run"] = float(bb1[key]["fail"].mean())
        ta5 = tails.analyze_tail(S5, tail_frac=max(1e-2, 1000 / len(S5)), censored_at=P["cens_at"])
        r["S5_gpd"] = ta5.get("gpd")
        r["S5_best_aic"] = ta5.get("best_aic")
        r["S5_pareto_alpha"] = ta5.get("fits", {}).get("pareto", {}).get("params", {}).get("alpha")
        out[str(key)] = r
        for S, lab in [(S1, "S_1: first solution"), (S5, "S_5: 5 solutions (stop_nconv=5)")]:
            x, sf = tails.survival(S)
            ok = sf > 0
            ax.loglog(x[ok], sf[ok], ".", ms=2, label=lab)
        ax.set_title(f"Relay-BP {P['label']} p={P['p']:g}", fontsize=8)
        ax.set_xlabel("BP iterations")
        ax.set_ylabel("P(S > x)")
        ax.legend(fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "exp2_relay_S1_vs_S5.png"))
    plt.close(fig)
    dump(out, "exp2_relay_nconv5.json")


def fig_param_trends(tail_rows):
    if not tail_rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for fam, mk in [("mwpm", "o-"), ("relay", "s-")]:
        rows = [r for r in tail_rows if r["family"] == fam]
        for p in sorted({r["p"] for r in rows}):
            rr = sorted([r for r in rows if r["p"] == p], key=lambda r: (r["d"], r["point"]))
            if fam == "mwpm":
                xs = [r["d"] for r in rr]
            else:
                xs = [BB_CODES[r["point"].strip("[]").replace(",", "_")].n for r in rr]
            lab = f"{'MWPM' if fam == 'mwpm' else 'Relay'} p={p:g}"
            axes[0].errorbar(xs, [r["xi"] for r in rr],
                             yerr=[[max(r["xi"] - r["xi_lo"], 0) for r in rr], [max(r["xi_hi"] - r["xi"], 0) for r in rr]],
                             fmt=mk, ms=3, capsize=2, label=lab)
            axes[1].plot(xs, [r.get("lognormal_sigma", np.nan) for r in rr], mk, ms=3, label=lab)
    axes[0].axhline(0, color="k", lw=0.7)
    axes[0].set_title(r"GPD shape $\xi$ above $u$ ($\xi>0$: Pareto-type)")
    axes[1].set_title(r"lognormal tail fit $\sigma$")
    for ax in axes:
        ax.set_xlabel("d (MWPM) / n (BB)")
        ax.set_xscale("log")
    axes[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "exp2_param_trends.png"))
    plt.close(fig)


def fig_buffers(buffers):
    if not buffers:
        return
    tags = list(buffers)
    fig, axes = plt.subplots(1, len(tags), figsize=(3.2 * len(tags), 3.1), squeeze=False, sharey=False)
    for ax, tag in zip(axes[0], tags):
        B = buffers[tag]
        best = B["best"]
        tbl = B["table"].get(best)
        if not isinstance(tbl, list):
            continue
        eps = [r["eps"] for r in tbl]
        for m, ls in [(1, "-"), (2, "--"), (3, ":")]:
            ax.loglog(eps, [r[f"med_m{m}"] for r in tbl], "k" + ls,
                      label=f"{'sliding (D/G/1)' if m == 1 else f'PWD m={m}'} [{best}]")
        for name in ("lognormal", "pareto"):
            t2 = B["table"].get(name)
            if isinstance(t2, list) and name != best:
                ax.loglog(eps, [r["med_m1"] for r in t2], color=COLORS[name], lw=0.8, label=f"D/G/1 [{name}]")
        ax.invert_xaxis()
        ax.set_xlabel(r"target $\varepsilon$")
        ax.set_title(B["title"], fontsize=7.5)
        ax.legend(fontsize=5.5)
    axes[0][0].set_ylabel("buffer (units of median decode work)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "exp2_buffer_vs_eps.png"))
    plt.close(fig)


def fig_lindley(lind):
    fig, axes = plt.subplots(1, len(lind), figsize=(3.6 * len(lind), 3.0), squeeze=False)
    for ax, (tag, res) in zip(axes[0], lind.items()):
        for rho, r in res.items():
            l, = ax.loglog(r["thresholds_med"], np.where(np.array(r["P_overflow"]) > 0, r["P_overflow"], np.nan),
                           "o-", ms=3, label=f"Lindley sim rho={rho}")
            ax.loglog(r["thresholds_med"], r["veraverbeke_empirical"], "--", color=l.get_color(), lw=1,
                      label=r"$\frac{\rho}{1-\rho}\bar F_I(x)$")
        ax.set_xlabel("buffer x (units of median S)")
        ax.set_ylabel("P(W > x)")
        ax.set_title(tag, fontsize=8)
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "exp2_lindley.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
