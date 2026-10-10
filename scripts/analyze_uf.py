"""Union-find (ldpc) work-distribution analysis on the surface-code data in data/uf.

Reuses the Experiment 2 machinery of scripts/analyze.py (tail fits, buffer sizing, diagnostic
figures) and compares each point with the MWPM data at the same (d, p).

Outputs: results/uf_*.json, results/uf_tail_fits.csv and figures/uf_*.png.
Usage: PYTHONPATH=. python3 scripts/analyze_uf.py
"""
import csv
import glob
import os
import re
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import relay_bp  # noqa: F401  (analyze.py imports it; not needed for the UF analysis)
except ImportError:
    import types
    sys.modules["relay_bp"] = types.ModuleType("relay_bp")
import analyze as an  # noqa: E402
from analyze import FIG, RES, MODEL_NAMES, dump, plt  # noqa: E402
from synbuf import exp1, tails  # noqa: E402
from synbuf.run_uf import COUNTERS  # noqa: E402

D_COLORS = {11: "#1f77b4", 15: "#d62728"}


def load_uf():
    pts = {}
    for f in sorted(glob.glob("data/uf/surface_d*_p*.npz")):
        m = re.search(r"d(\d+)_p([0-9.e-]+)\.npz", f)
        d, p = int(m.group(1)), float(m.group(2))
        z = np.load(f)
        pts[(d, p)] = dict(S=z["S"], fail=z["fail"], cnt={k: z[k] for k in COUNTERS},
                           label=f"surface d={d}", d=d, p=p, decoder="UF")
    return pts


def load_mwpm(keys):
    out = {}
    for d, p in keys:
        f = f"data/mwpm/surface_d{d}_p{p:g}.npz"
        if os.path.exists(f):
            z = np.load(f)
            out[(d, p)] = dict(S=z["S"], fail=z["fail"])
    return out


def point_summary(P, M):
    S = P["S"].astype(np.float64)
    k = int(P["fail"].sum())
    n = len(S)
    row = dict(d=P["d"], p=P["p"], shots=n, ler=k / n, ler_ci=exp1.wilson(k, n), fails=k,
               mean=S.mean(), median=float(np.median(S)), std=S.std(), max=S.max(),
               max_over_median=S.max() / np.median(S),
               q={f"1e-{j}": float(tails.quantile_tail(S, 10.0 ** -j)) for j in range(1, 7) if 10.0 ** -j * n >= 1},
               counters_mean={c: float(P["cnt"][c].mean()) for c in COUNTERS},
               counters_max={c: int(P["cnt"][c].max()) for c in COUNTERS})
    # robustness of the work metric: rank correlation with finer-grained cost models
    fine = (P["cnt"]["bits_added"].astype(np.float64) + P["cnt"]["merge_moves"])
    row["spearman_S_vs_bits_plus_moves"] = float(stats.spearmanr(S, fine).statistic)
    row["spearman_S_vs_defects"] = float(stats.spearmanr(S, P["cnt"]["clusters_created"]).statistic)
    # node-level cost W = bits added + nodes copied by merges (ldpc merges by copying sets)
    ta = tails.analyze_tail(fine, tail_frac=1e-3)
    row["node_cost"] = dict(mean=fine.mean(), median=float(np.median(fine)), max=fine.max(),
                            max_over_median=fine.max() / np.median(fine), best_aic=ta.get("best_aic"),
                            xi=(ta.get("gpd") or {}).get("xi"), xi_ci=(ta.get("gpd") or {}).get("xi_ci"))
    if M is not None:
        Sm = M["S"].astype(np.float64)
        km = int(M["fail"].sum())
        row["mwpm"] = dict(ler=km / len(Sm), ler_ci=exp1.wilson(km, len(Sm)), fails=km, mean=Sm.mean(),
                           median=float(np.median(Sm)), max=Sm.max(), max_over_median=Sm.max() / np.median(Sm))
    return row


def fig_ccdf(uf, mw):
    ps = sorted({k[1] for k in uf})
    fig, axes = plt.subplots(1, len(ps), figsize=(4.2 * len(ps), 3.4), sharey=True, squeeze=False)
    for ax, p in zip(axes[0], ps):
        for d in sorted({k[0] for k in uf if k[1] == p}):
            for dec, src, ls in [("UF", uf, "-"), ("MWPM", mw, "--")]:
                if (d, p) not in src:
                    continue
                S = src[(d, p)]["S"].astype(np.float64)
                x, sf = tails.survival(S)
                ok = sf > 0
                ax.loglog(x[ok] / np.median(S), sf[ok], ls, color=D_COLORS.get(d, "k"), lw=1.2,
                          label=f"{dec} d={d}")
        ax.set_xlabel("S / median(S)")
        ax.set_title(f"p = {p:g}")
        ax.legend(fontsize=7)
    axes[0][0].set_ylabel(r"$P(S > x)$")
    fig.suptitle("Work distribution tails: union-find (ldpc) vs MWPM, normalized by median", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "uf_ccdf_vs_mwpm.png"))
    plt.close(fig)


def main():
    uf = load_uf()
    if not uf:
        sys.exit("no data in data/uf")
    mw = load_mwpm(uf)
    print("UF points:", sorted(uf), " MWPM matches:", sorted(mw))

    summary = {f"d{d}_p{p:g}": point_summary(P, mw.get((d, p))) for (d, p), P in sorted(uf.items())}

    # ---------------- tails (same settings as MWPM in analyze.py)
    tail_rows, buffers = [], {}
    for (d, p), P in sorted(uf.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        frac = 1e-3 if len(P["S"]) * 1e-3 >= 2000 else 2000 / len(P["S"])
        ta = an.tail_analysis(P, frac)
        ta2 = tails.analyze_tail(P["S"], tail_frac=frac / 10)
        row = dict(decoder="UF", d=d, p=p, n=ta["n"], mean=ta["mean"], median=ta["median"], u=ta["u"],
                   n_exceed=ta["n_exceed"], best_aic=ta.get("best_aic"),
                   xi=ta.get("gpd", {}).get("xi"), xi_lo=ta.get("gpd", {}).get("xi_ci", [None, None])[0],
                   xi_hi=ta.get("gpd", {}).get("xi_ci", [None, None])[1],
                   xi_high_u=(ta2.get("gpd") or {}).get("xi"), best_aic_high_u=ta2.get("best_aic"),
                   n_exceed_high_u=ta2["n_exceed"])
        for name in MODEL_NAMES:
            if name in ta["fits"]:
                for pk, pv in ta["fits"][name]["params"].items():
                    row[f"{name}_{pk}"] = pv
                row[f"{name}_aic"] = ta["aic"][name]
        for pair, (R, z, pv) in ta["vuong"].items():
            row[f"vuong_{pair}_z"] = z
            row[f"vuong_{pair}_p"] = pv
        tail_rows.append(row)
        print(f"[tail] UF d={d} p={p:g}: mean={row['mean']:.1f} best={row['best_aic']} xi={row['xi']:.3f} "
              f"(high-u xi={row['xi_high_u']})", flush=True)
        tag = f"{d}_{p:g}"
        title = f"Union-find (ldpc) surface d={d} p={p:g}"
        an.fig_diagnostics(P, ta, f"uf_diag_{tag}.png", title)
        an.fig_extrapolation(P, ta, f"uf_extrap_{tag}.png", title)
        buffers[f"d{d}_p{p:g}"] = dict(title=title, mean=float(P["S"].mean()), median=float(np.median(P["S"])),
                                       best=ta.get("best_aic"), table=an.buffer_table(P, ta))
    dump(tail_rows, "uf_tail_fits.json")
    cols = sorted({k for r in tail_rows for k in r}, key=lambda c: (c not in ("decoder", "d", "p"), c))
    with open(os.path.join(RES, "uf_tail_fits.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(tail_rows)
    dump(buffers, "uf_buffers.json")
    dump(summary, "uf_summary.json")
    fig_ccdf(uf, mw)

    for k, s in summary.items():
        m = s.get("mwpm", {})
        print(f"{k}: LER={s['ler']:.2e} (MWPM {m.get('ler', float('nan')):.2e})  S mean={s['mean']:.1f} "
              f"median={s['median']:.0f} max={s['max']:.0f} max/med={s['max_over_median']:.2f} "
              f"(MWPM max/med={m.get('max_over_median', float('nan')):.2f})  node-cost max/med="
              f"{s['node_cost']['max_over_median']:.1f} xi={s['node_cost']['xi']:.3f}")


if __name__ == "__main__":
    main()
