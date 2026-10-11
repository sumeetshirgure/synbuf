"""Tail class of decoder work: union-find (ldpc, peeling) vs MWPM (Sparse Blossom), surface code.

Data: data/uf (synbuf.run_uf) and the matching points of data/mwpm (synbuf.run_mwpm); the MWPM
files must carry the full counter set (any run with the current PyMatching patch).

Every counter counts bounded-cost operations, so decode time is within constant factors of the
counter sum S, and the tail class of S is that of its heaviest-tailed component whatever the
per-counter costs. The analysis therefore fits the tail of S and of each component group (dual /
primal / extraction for MWPM; growth / merge / validity / peel for UF) with a GPD and reports the
shape xi (profile-likelihood 95% CI, LR test of xi = 0); no per-operation cost constants are used.

Outputs: results/uf_summary.json, results/uf_tail_fits.{json,csv}, figures/uf_*.png.
Usage: PYTHONPATH=. python3 scripts/analyze_uf.py
"""
import csv
import os
import sys

import matplotlib.ticker
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import relay_bp  # noqa: F401  (analyze.py imports it; not needed here)
except ImportError:
    import types
    sys.modules["relay_bp"] = types.ModuleType("relay_bp")
import analyze as an  # noqa: E402
from analyze import FIG, RES, dump, plt  # noqa: E402
from synbuf import exp1, tails  # noqa: E402
from synbuf.run_mwpm import COUNTERS as MWPM_COUNTERS  # noqa: E402
from synbuf.run_uf import COUNTERS as UF_COUNTERS  # noqa: E402

GROUPS = {
    "MWPM": {
        "dual": ["pq_pops", "queue_moves", "flood_events", "node_reschedules", "neighbor_scans"],
        "primal": ["mwpm_events", "tree_ops", "region_ops", "node_visits"],
        "extraction": ["shatter_ops", "extract_ops"],
    },
    "UF": {
        "growth": ["clusters_created", "growth_rounds", "cluster_grows", "growth_row_visits", "bit_add_calls",
                   "growth_col_visits"],
        "merge": ["merges", "merge_moves"],
        "validity": ["cluster_scans", "sort_compares"],
        "peel": ["peel_visits", "peel_steps", "cluster_solves"],
    },
}
COUNTERS = {"MWPM": MWPM_COUNTERS, "UF": UF_COUNTERS}
DIRS = {"MWPM": "data/mwpm", "UF": "data/uf"}
POINTS = [(d, p) for d in (11, 15) for p in (0.001, 0.003)]
D_COLORS = {11: "#1f77b4", 15: "#d62728"}
GROUP_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
TAIL_FRAC = 1e-3

for dec, groups in GROUPS.items():
    assert sorted(c for g in groups.values() for c in g) == sorted(COUNTERS[dec]), dec


def load(dec):
    pts = {}
    for d, p in POINTS:
        f = os.path.join(DIRS[dec], f"surface_d{d}_p{p:g}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f)
        missing = [c for c in COUNTERS[dec] if c not in z]
        if missing:
            sys.exit(f"{f} lacks counters {missing}; rerun it with the current decoder patch")
        P = dict(S=z["S"].astype(np.float64), fail=z["fail"], d=d, p=p, label=f"surface d={d}")
        for g, cs in GROUPS[dec].items():
            P[g] = sum(z[c].astype(np.float64) for c in cs)
        P["counters"] = {c: z[c] for c in COUNTERS[dec]}
        pts[(d, p)] = P
    return pts


def fit(x, frac=TAIL_FRAC):
    """GPD tail fit of samples x at the (1-frac) quantile, and again at a 10x higher threshold."""
    fr = frac if len(x) * frac >= 2000 else 2000 / len(x)
    ta = tails.analyze_tail(x, tail_frac=fr)
    hi = tails.analyze_tail(x, tail_frac=fr / 10)
    g, g2 = ta["gpd"], hi["gpd"]
    return dict(u=ta["u"], n_exceed=ta["n_exceed"], xi=g["xi"], xi_ci=g["xi_ci"], scale=g["scale"],
                lr_exp=g["lr_exp"], u_high=hi["u"], n_exceed_high=hi["n_exceed"], xi_high=g2["xi"],
                xi_ci_high=g2["xi_ci"], lr_exp_high=g2["lr_exp"])


def describe(x):
    med = float(np.median(x))
    return dict(mean=float(x.mean()), median=med, std=float(x.std()), max=float(x.max()),
                max_over_median=float(x.max() / med) if med > 0 else None,
                q={f"1e-{j}": float(tails.quantile_tail(x, 10.0 ** -j)) for j in range(1, 7) if 10.0 ** -j * len(x) >= 1})


def fig_ccdf(data):
    ps = sorted({k[1] for P in data.values() for k in P})
    fig, axes = plt.subplots(1, len(ps), figsize=(4.2 * len(ps), 3.4), sharey=True, squeeze=False)
    for ax, p in zip(axes[0], ps):
        for d in sorted({k[0] for P in data.values() for k in P if k[1] == p}):
            for dec, ls in [("UF", "-"), ("MWPM", "--")]:
                if (d, p) not in data[dec]:
                    continue
                S = data[dec][(d, p)]["S"]
                x, sf = tails.survival(S)
                ok = sf > 0
                ax.loglog(x[ok] / np.median(S), sf[ok], ls, color=D_COLORS.get(d, "k"), lw=1.2, label=f"{dec} d={d}")
        ax.set_xlabel("S / median(S)")
        ax.set_title(f"p = {p:g}")
        ax.legend(fontsize=7)
    axes[0][0].set_ylabel(r"$P(S > x)$")
    fig.suptitle("Total work S (sum of all bounded-cost counters), normalized by median", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "uf_ccdf_vs_mwpm.png"))
    plt.close(fig)


def fig_components(dec, pts):
    keys = sorted(pts, key=lambda k: (k[1], k[0]))
    fig, axes = plt.subplots(1, len(keys), figsize=(3.3 * len(keys), 3.2), sharey=True, squeeze=False)
    for ax, key in zip(axes[0], keys):
        P = pts[key]
        for (g, color) in list(zip(GROUPS[dec], GROUP_COLORS)) + [("S", "k")]:
            x = P[g]
            if x.max() <= 0:
                continue
            xs, sf = tails.survival(x)
            ok = (sf > 0) & (xs > 0)
            ax.loglog(xs[ok] / np.median(x[x > 0]), sf[ok], color=color, lw=1.6 if g == "S" else 1.0,
                      ls="-" if g == "S" else "--", label="S (total)" if g == "S" else g)
        ax.set_title(f"d={key[0]}, p={key[1]:g}", fontsize=8)
        ax.set_xlabel("x / median")
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axes[0][0].set_ylabel(r"$P(\cdot > x)$")
    axes[0][0].legend(fontsize=6)
    fig.suptitle(f"{dec}: work components, each normalized by its median", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"uf_components_{dec.lower()}.png"))
    plt.close(fig)


def main():
    data = {dec: load(dec) for dec in GROUPS}
    for dec, pts in data.items():
        print(dec, "points:", sorted(pts))
        if not pts:
            sys.exit(f"no data in {DIRS[dec]}")

    summary, rows = {}, []
    for dec, pts in data.items():
        for (d, p), P in sorted(pts.items(), key=lambda kv: (kv[0][1], kv[0][0])):
            n = len(P["S"])
            k = int(P["fail"].sum())
            entry = dict(decoder=dec, d=d, p=p, shots=n, ler=k / n, ler_ci=exp1.wilson(k, n), fails=k,
                         S=describe(P["S"]),
                         components={g: describe(P[g]) for g in GROUPS[dec]},
                         share_of_mean={g: float(P[g].mean() / P["S"].mean()) for g in GROUPS[dec]},
                         counters_mean={c: float(v.mean()) for c, v in P["counters"].items()},
                         counters_max={c: int(v.max()) for c, v in P["counters"].items()})
            fits = {}
            for metric in ["S"] + list(GROUPS[dec]):
                x = P[metric]
                if x.max() <= 0:
                    continue
                fr = fit(x)
                fits[metric] = fr
                desc = entry["S"] if metric == "S" else entry["components"][metric]
                rows.append(dict(decoder=dec, d=d, p=p, metric=metric, median=desc["median"], mean=desc["mean"],
                                 max_over_median=desc["max_over_median"], u=fr["u"], n_exceed=fr["n_exceed"],
                                 xi=fr["xi"], xi_lo=fr["xi_ci"][0], xi_hi=fr["xi_ci"][1], scale=fr["scale"],
                                 lr_exp_stat=fr["lr_exp"]["stat"], lr_exp_p=fr["lr_exp"]["p"],
                                 u_high=fr["u_high"], n_exceed_high=fr["n_exceed_high"], xi_high=fr["xi_high"],
                                 xi_lo_high=fr["xi_ci_high"][0], xi_hi_high=fr["xi_ci_high"][1],
                                 lr_exp_p_high=fr["lr_exp_high"]["p"]))
                if metric == "S":
                    an.fig_ccdf(x, f"uf_ccdf_{dec.lower()}_{d}_{p:g}.png",
                                f"{dec} surface d={d} p={p:g}: total work S", xlabel="S (operations)")
            comp = [g for g in GROUPS[dec] if g in fits]
            heaviest = max(comp, key=lambda g: fits[g]["xi"])
            entry["tail"] = dict(S=fits["S"], heaviest_component=heaviest, components={g: fits[g] for g in comp})
            summary[f"{dec}_d{d}_p{p:g}"] = entry
            cx = ", ".join(f"{g} {fits[g]['xi']:+.3f}" for g in comp)
            fS = fits["S"]
            print(f"{dec:4s} d={d} p={p:g}: LER={entry['ler']:.2e} S median={entry['S']['median']:.0f} "
                  f"max/med={entry['S']['max_over_median']:.2f} xi(S)={fS['xi']:+.3f} "
                  f"[{fS['xi_ci'][0]:+.3f},{fS['xi_ci'][1]:+.3f}] LR(xi=0) p={fS['lr_exp']['p']:.2g} | {cx}",
                  flush=True)
        fig_components(dec, pts)
    fig_ccdf(data)

    dump(summary, "uf_summary.json")
    dump(rows, "uf_tail_fits.json")
    cols = list(rows[0])
    with open(os.path.join(RES, "uf_tail_fits.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
