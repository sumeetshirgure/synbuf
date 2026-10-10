# synbuf: Results Report (Experiments 1 and 2)

*Status as of 2026-10-10. Describes the runs actually performed and the results currently on disk in `data/`, `results/` and `figures/`. The plan these runs implement is in [EXPERIMENTS.md](EXPERIMENTS.md); deviations from it are listed in [§2](#2-what-was-and-was-not-run). The earlier MWPM runs ($d=3,7,11,15$, legacy work metric) and every result built on them were removed on 2026-10-10. MWPM is now covered only in [§7](#7-union-find-decoder-ldpc-vs-mwpm-full-operation-counts), with full operation counts; §4–5 cover Relay-BP only.*

## 1. Summary of findings

1. **Capping Relay-BP blows up the logical error rate; it does not produce a flat error floor.** With a fixed work budget $c$, the capped LER $P_L(c)$ tracks the uncapped LER at low $p$ and then rises far faster than it as $p$ grows. For [[144,12,12]] with $c=200$ iterations, $P_L(c)=2\times10^{-6}$, $5.5\times10^{-5}$ and $1.5\times10^{-2}$ at $p=5\times10^{-4}$, $10^{-3}$ and $3\times10^{-3}$, while the uncapped LER is 0, 0 and $5.3\times10^{-4}$. The excess is set by $P(S>c)$, which is strongly $p$-dependent. So the cap acts as an additional, rapidly growing failure mechanism rather than as a constant offset.
2. **Relay-BP work has a much heavier, power-law tail that grows heavier with $p$.** $\xi$ ranges from about 0.05–0.4 at $p=5\times10^{-4}$, through 0.34–0.94 at $10^{-3}$, to above 1 at $3\times10^{-3}$ for the three codes with $d\ge10$ ($\alpha\approx0.74$–$1.01$). A fit at a 10× higher threshold gives $\xi\approx0.84$–$1.15$ there, so the tail sits at or near the infinite-mean boundary rather than clearly beyond it. With $\alpha\le1$ the leading-order buffer formulas diverge, and a D/G/1 Lindley simulation of [[144,12,12]] at $p=3\times10^{-3}$ spends about a quarter of its time draining single huge decodes even at load $\rho=0.5$. Spare servers (parallel windows) shrink the Pareto-tail buffer by orders of magnitude where the formulas still apply.
3. **With all bounded-cost operations counted, union-find has a lighter tail than MWPM, and MWPM's heaviest tail comes from its primal (matching) side** ([§7](#7-union-find-decoder-ldpc-vs-mwpm-full-operation-counts)). At $d\in\{11,15\}$ and $p\in\{10^{-3},3\times10^{-3}\}$, MWPM's total work has $\xi\approx0.17$–$0.31$ (Pareto-type) and a maximum of 8–12× the median. Its primal work is only about 10% of the mean but has the heaviest tail ($\xi\approx0.21$–$0.45$, rising with $p$ and $d$). Peeling UF has $\xi\approx0.14$ at $p=10^{-3}$, $\xi\approx0$ at $p=3\times10^{-3}$, and a maximum of 4–6× the median. This UF variant is much less accurate than MWPM: about 36× and 210× higher LER at $p=3\times10^{-3}$ for $d=11$ and $d=15$.

## 2. What was and was not run

### Run

| Item | What was done |
|---|---|
| Relay-BP data | BB codes [[72,12,6]], [[90,8,10]], [[108,8,10]], [[144,12,12]], $p\in\{5\times10^{-4},10^{-3},3\times10^{-3}\}$, $10^6$ shots per point (12 points) |
| MWPM and UF data | Rotated surface code memory, $d\in\{11,15\}$, $p\in\{10^{-3},3\times10^{-3}\}$, $10^6$ shots per point for each decoder (§7). The earlier MWPM runs ($d\in\{3,7,11,15\}$ at three values of $p$, legacy metric) were removed. |
| Instrumented decoders | PyMatching, Relay-BP and ldpc (union-find), built from pinned commits with the patches in `third_party/` (`scripts/build_decoders.sh`) |
| Relay-BP parameters | Per-code tuned values from `results/relay_params.json` (tuned earlier; see `results/relay_tuning_notes.md`) |
| Experiment 1 (Relay-BP) | Offline $P_L(c)$ for an absolute cap ladder, the deployed Relay-BP limit, and calibrated caps $c_3,c_4,c_5$ |
| Experiment 2 (Relay-BP) | Peaks-over-threshold tail fits (GPD, exponential, lognormal, Pareto), Vuong tests, AIC ranking and censored likelihood for every point; buffer tables for 5 selected points; an empirical log-log survival plot for every point; D/G/1 Lindley simulation at 1 point |

### Not run, or run differently from EXPERIMENTS.md

| Planned in EXPERIMENTS.md | Status |
|---|---|
| MWPM: odd $d$ from 3 to about 15, three values of $p$ | Only $d=11,15$ at $p\in\{10^{-3},3\times10^{-3}\}$ (§7); no MWPM capping experiment |
| Uniform depolarizing **or** SI1000 noise | Uniform circuit-level depolarizing only |
| Validate BB-code LERs against published values | Not part of this run |
| Optional coprime BB codes | Not included. [[288,12,18]] has tuned parameters but no data in this run |
| Retune Relay-BP per code | Not re-run; the existing tuned parameters were reused (the `tune` stage skips when `results/relay_params.json` exists). Tuning used $p=0.004$–$0.005$ and a 300-leg limit |
| Report time to $S>1$ solutions separately (`nconv5` stage, `stop_nconv=5`) | Not run; only time to first solution ($S_1$) |
| Raise the leg limit until censoring is far below the fitted tail probabilities | Not done. Up to $1.3\times10^{-4}$ of shots are censored at $p=3\times10^{-3}$; the limit (600,080 iterations) is already far beyond any practical cap, so censored shots enter the fits through a censored likelihood instead |
| Caps defined as a hardware budget (ops per µs at an assumed clock) | Caps are absolute operation counts on a 1-2-5 grid, not tied to a specific clock or hardware budget |
| Headline figure: LER vs $d$ at fixed $c$ | Replaced by $P_L(c)$ vs $p$ at fixed $(d,c)$ (§4.3) |
| $P_L(c)$ vs $c$, and $P(F_\infty\mid S>c)$ vs $c$ | Removed from the analysis |
| Tail-parameter trends with $d$ and $p$ (figure) | Parameters are reported in tables (§5.3); the trend figure was removed |
| Buffer size vs $\varepsilon$ (figure) | Buffer values are in `results/exp2_buffers.json` and summarized in §5.4; the figure was removed |
| PWD (D/G/k) buffer | Leading-order formulas only; no multi-server simulation |
| Lindley check | D/G/1 only, at 1 point; numbers in `results/exp2_lindley.json`, figure removed |
| Diagnostic plots and fitted-tail extrapolation plots | Removed; each point has an empirical log-log survival plot without fitted curves |

## 3. Common setup

**Circuits.** Generated with Stim (`synbuf/circuits.py`).
- Surface code: Stim's rotated memory circuit with $d$ rounds, with `after_clifford_depolarization`, `before_round_data_depolarization` and `before_measure_flip_probability` all set to $p$.
- BB codes: the Bravyi et al. (Nature 627, 778, 2024) syndrome cycle with uniform circuit-level noise, $d$ rounds.

**Window.** One $d$-round memory block is decoded as one window. All LERs are per window, and for BB codes per window of the whole block (any logical error among the $k$ logical qubits counts as one failure).

**Work metric $S$.**
- MWPM and union-find: the sum of all bounded-cost operation counters (§7).
- Relay-BP: $S$ = total BP iterations up to and including the first converged leg (time to first solution, `stop_nconv=1`). Each decode runs 80 pre-iterations, then up to 10,000 legs of up to 60 iterations each, so the leg limit is $80+10{,}000\times60=600{,}080$ iterations. A shot that reaches it without converging is **censored** and counted as a logical failure.

**Sampling.** Shots are split into chunks with per-chunk seeds and decoded by 10 worker processes, so results do not depend on the worker count. Data: `data/relay/bb{code}_p{p}/`, `data/mwpm/surface_d{d}_p{p}.npz` and `data/uf/surface_d{d}_p{p}.npz`.

## 4. Experiment 1: decoder caps and the LER (Relay-BP)

### 4.1 Method

Every shot is decoded to completion and records $(S, \mathbb 1[F_\infty])$. For any cap $c$ the capped LER is computed offline, with an aborted decode counted as a failure:

$$P_L(c) = P(F_\infty,\ S\le c) + P(S>c).$$

### 4.2 Caps

- **Absolute ladder.** Caps on a 1-2-5 grid from 2 to 500,000 iterations (from the smallest median $S$ to the largest observed $S$, kept below the leg limit). Each cap is held fixed while the code and $p$ vary.
- **Deployed Relay-BP limit.** $80 + 300\times60 = 18{,}080$ iterations (300 legs, the limit used during tuning).
- **Calibrated caps.** $c_j$ with $P(S>c_j)=10^{-j}$ at the reference point ([[144,12,12]], $p=3\times10^{-3}$): $c_3=2{,}806$ and $c_4=66{,}456$. $c_5$ cannot be calibrated, because the censored fraction there ($2.6\times10^{-5}$) exceeds $10^{-5}$. These caps are listed in `results/exp1_relay.{json,csv}` but not plotted.

### 4.3 Results

Figure: [`figures/exp1_ler_vs_p_relay.png`](figures/exp1_ler_vs_p_relay.png). Each panel is one code. It shows:
- the uncapped LER (dashed black);
- up to five ladder caps that bind somewhere on the $p$ grid but cap fewer than 5% of shots at the lowest $p$;
- the deployed limit (purple).

**Relay-BP, $P_L(c)$ per window** (deployed limit = 18,080 iterations; a value of 0 means no failures in $10^6$ shots):

| Code | $p$ | median $S$ | uncapped | $c=200$ | $c=2000$ | deployed | $c=20{,}000$ | $c=200{,}000$ |
|---|---|---|---|---|---|---|---|---|
| [[72,12,6]] | 5e-4 | 1 | 2.7e-5 | 2.7e-5 | 2.7e-5 | 2.7e-5 | 2.7e-5 | 2.7e-5 |
| [[72,12,6]] | 1e-3 | 2 | 2.89e-4 | 3.01e-4 | 2.89e-4 | 2.89e-4 | 2.89e-4 | 2.89e-4 |
| [[72,12,6]] | 3e-3 | 8 | 1.73e-2 | 1.93e-2 | 1.73e-2 | 1.73e-2 | 1.73e-2 | 1.73e-2 |
| [[90,8,10]] | 5e-4 | 2 | 1.0e-6 | 2.6e-5 | 4.0e-6 | 3.0e-6 | 3.0e-6 | 2.0e-6 |
| [[90,8,10]] | 1e-3 | 4 | 4.0e-6 | 2.54e-4 | 5.1e-5 | 2.4e-5 | 2.4e-5 | 1.1e-5 |
| [[90,8,10]] | 3e-3 | 11 | 2.43e-3 | 2.30e-2 | 5.20e-3 | 3.05e-3 | 3.01e-3 | 2.54e-3 |
| [[108,8,10]] | 5e-4 | 2 | 0 | 2.7e-5 | 3.0e-6 | 3.0e-6 | 3.0e-6 | 0 |
| [[108,8,10]] | 1e-3 | 4 | 4.0e-6 | 2.50e-4 | 2.8e-5 | 1.3e-5 | 1.3e-5 | 7.0e-6 |
| [[108,8,10]] | 3e-3 | 10 | 1.23e-3 | 2.31e-2 | 3.76e-3 | 1.82e-3 | 1.79e-3 | 1.32e-3 |
| [[144,12,12]] | 5e-4 | 4 | 0 | 2.0e-6 | 0 | 0 | 0 | 0 |
| [[144,12,12]] | 1e-3 | 6 | 0 | 5.5e-5 | 3.0e-6 | 1.0e-6 | 1.0e-6 | 1.0e-6 |
| [[144,12,12]] | 3e-3 | 15 | 5.33e-4 | 1.54e-2 | 1.73e-3 | 7.32e-4 | 7.19e-4 | 5.53e-4 |

The uncapped Relay-BP LER counts censored shots (those reaching 600,080 iterations) as failures.

### 4.4 Findings

- **Blow-up, not a floor.** For a fixed cap, $P_L(c)$ converges to the uncapped LER as $p$ falls and rises much more steeply than it as $p$ grows. For [[144,12,12]] at $c=200$: $2\times10^{-6}$, $5.5\times10^{-5}$ and $1.5\times10^{-2}$ at $p=5\times10^{-4}$, $10^{-3}$ and $3\times10^{-3}$ (29× the uncapped LER at the highest $p$). A floor would require $P(S>c)$ to be insensitive to $p$, but it is not: the cap opens a new failure channel whose rate grows quickly with $p$.
- **Relay-BP needs caps far above its median.** Because of its heavy tail, caps of $10^2$–$10^4\times$ the median still add measurable error at $p\ge10^{-3}$ for the $d\ge10$ codes. The deployed 300-leg limit raises the LER by 1.3–1.5× at $p=3\times10^{-3}$ and by about 3–6× at $p=10^{-3}$ for [[90,8,10]] and [[108,8,10]], though the low-$p$ counts are small (4 uncapped vs 24 capped failures for [[90,8,10]]). [[72,12,6]] is insensitive to every cap of 2,000 iterations or more because its LER is high.
- **The distance trend is not established.** The earlier MWPM runs suggested that a fixed cap hurts more at larger $d$. Those runs have been removed, and the BB codes do not form a distance family, so this report makes no claim about how capping scales with $d$.
- **Statistics.** At $p\le10^{-3}$ for the $d\ge10$ BB codes, the uncapped LER is at or below the $10^{-6}$ resolution of $10^6$ shots, so ratios at those points are lower bounds or rest on a handful of events. Error bars in the figures are 1σ Wilson intervals.

## 5. Experiment 2: tail class of the work distribution (Relay-BP)

### 5.1 Method

Implemented in `synbuf/tails.py` and `scripts/analyze.py`, on the uncapped Experiment 1 data.
- **Threshold.** $u$ is the empirical $(1-q)$ quantile of $S$ with $q=10^{-2}$ (about 10,000 exceedances). A second fit at $q/10$ checks threshold stability. (§7 uses $q=10^{-3}$.)
- **Fits.** Left-truncated MLE fits above $u$ of exponential, lognormal and Pareto tails, plus a generalized Pareto (GPD) fit with a bootstrap 95% interval on its shape $\xi$. $\xi>0$ indicates a Pareto-type tail; $\xi\approx0$ the Gumbel domain (exponential or lognormal).
- **Model selection.** AIC ranking and pairwise Vuong likelihood-ratio tests.
- **Censoring.** Relay-BP shots at the leg limit enter the likelihood as right-censored observations ($S>600{,}080$).
- **Buffers.** The fitted tail is spliced onto the empirical bulk below $u$. For sliding windows (D/G/1), solve $\frac{\rho}{1-\rho}\bar F_I(x)=\varepsilon$ (Veraverbeke), where $\bar F_I$ is the integrated tail. For parallel windows with $m$ spare servers, solve $\bar F_I(x)^m=\varepsilon$ (Foss–Korshunov, leading order). $\rho=0.5$.
- **Lindley check.** A D/G/1 queue simulated for $2\times10^8$ windows by resampling the empirical $S$, at $\rho=0.5$ and $0.8$.

### 5.2 Figures

- Empirical survival function $P(S>x)$ on log-log axes for every point (12 plots, no fitted curves): `figures/exp2_ccdf_relay_*.png`.

### 5.3 Tail fits

From `results/exp2_tail_fits.csv`. "Censored" is the number of shots at the leg limit.

| Code | $p$ | $u$ | Best (AIC) | $\xi$ [95% CI] | Pareto $\alpha$ | Censored |
|---|---|---|---|---|---|---|
| [[72,12,6]] | 5e-4 | 11 | lognormal | 0.05 [0.01, 0.08] | 3.0 | 0 |
| [[72,12,6]] | 1e-3 | 17 | Pareto | 0.34 [0.32, 0.37] | 3.2 | 0 |
| [[72,12,6]] | 3e-3 | 111 | lognormal | 0.53 [0.50, 0.56] | 1.6 | 0 |
| [[90,8,10]] | 5e-4 | 15 | Pareto | 0.23 [0.19, 0.28] | 3.6 | 0 |
| [[90,8,10]] | 1e-3 | 21 | Pareto | 0.80 [0.77, 0.85] | 2.4 | 0 |
| [[90,8,10]] | 3e-3 | 539 | Pareto | 1.42 [1.38, 1.46] | 0.74 | 117 |
| [[108,8,10]] | 5e-4 | 11 | Pareto | 0.42 [0.38, 0.46] | 3.4 | 0 |
| [[108,8,10]] | 1e-3 | 15 | Pareto | 0.94 [0.90, 0.98] | 2.2 | 1 |
| [[108,8,10]] | 3e-3 | 461 | Pareto | 1.46 [1.42, 1.50] | 0.76 | 128 |
| [[144,12,12]] | 5e-4 | 18 | lognormal | 0.09 [0.05, 0.11] | 4.2 | 0 |
| [[144,12,12]] | 1e-3 | 24 | Pareto | 0.48 [0.45, 0.51] | 3.3 | 0 |
| [[144,12,12]] | 3e-3 | 293 | Pareto | 1.07 [1.03, 1.11] | 1.01 | 26 |

### 5.4 Buffer estimates at $\varepsilon=10^{-12}$

From `results/exp2_buffers.json`, in units of the median decode work, at $\rho=0.5$, using the AIC-best model. "Diverges" means the integrated tail is infinite ($\alpha\le1$ or $\xi\ge1$), or the solver reached its $10^{40}$ upper bound.

| Point | Model | Sliding (D/G/1) | PWD, $m=2$ | PWD, $m=3$ | GPD, D/G/1 |
|---|---|---|---|---|---|
| Relay [[144,12,12]], $p=10^{-3}$ | Pareto | $9.5\times10^4$ | 247 | 34 | $4\times10^9$ |
| Relay [[72,12,6]], $p=3\times10^{-3}$ | lognormal | $5.1\times10^7$ | $1.2\times10^5$ | $8.0\times10^3$ | $5.3\times10^{13}$ |
| Relay [[90,8,10]], [[108,8,10]], [[144,12,12]], $p=3\times10^{-3}$ | Pareto | diverges | diverges | diverges | diverges |

### 5.5 Lindley check (D/G/1)

From `results/exp2_lindley.json`. $P(W>x)$ is the fraction of windows whose backlog exceeds $x$, with $x$ in multiples of the median $S$.

| Point | $\rho$ | $x=5$ | $x=50$ | $x=500$ | $x=5000$ |
|---|---|---|---|---|---|
| Relay [[144,12,12]], $p=3\times10^{-3}$ | 0.5 | 0.40 | 0.38 | 0.34 | 0.25 |
| Relay [[144,12,12]], $p=3\times10^{-3}$ | 0.8 | 0.74 | 0.71 | 0.68 | 0.59 |

For Relay-BP, each decode that reaches the leg limit (about 1 in 40,000 windows) takes about 8,000 windows to drain, so the queue spends about a quarter of its time recovering from one such decode. This is the single-big-jump regime in its extreme form. The Veraverbeke asymptotic exceeds 1 at these thresholds and is not valid here.

### 5.6 Findings

- **Relay-BP.**
  - The tail grows heavier with $p$. At every $p\ge10^{-3}$ it is heavier than MWPM's tail at the surface-code points of §7 (MWPM total work: $\xi\approx0.17$–$0.31$).
  - At $p=3\times10^{-3}$ the three $d\ge10$ codes have $\alpha\le1$ over the observed range. The single-server buffer has no finite leading-order size, and the simulated queue is dominated by individual huge decodes.
  - Where $\alpha>1$ ([[144,12,12]] at $p=10^{-3}$), spare servers cut the buffer from about $10^5$ to $34$ median decodes ($m=3$). This matches the prediction that parallel windows help most for a fat tail.
- **Threshold caveat.** For Relay-BP the threshold $u$ sits only a few times above the median at $p\le10^{-3}$ (11–24 iterations), so those fits describe the top 1% of the distribution rather than only the far tail. The higher-threshold fits in `exp2_tail_fits.csv` (`xi_high_u`, `best_aic_high_u`) also give $\xi>0$ at every point but move $\xi$ by up to about 1 at some points (for example [[108,8,10]] at $p=5\times10^{-4}$: 0.42 vs 1.46). At $p=3\times10^{-3}$ the higher threshold gives $\xi=0.86$, $0.84$ and $1.15$ for [[90,8,10]], [[108,8,10]] and [[144,12,12]], and prefers lognormal for the first two. So $\xi>1$ at that $p$ is threshold-sensitive; the robust statement is that the tail is at or near the infinite-mean boundary.
- **Censoring caveat.** Relay-BP censoring is at most $1.3\times10^{-4}$ of shots. It affects only the tail beyond 600,080 iterations, which is irrelevant for any practical cap but means buffer estimates that extrapolate past that point are effectively bounded by it.

## 6. Reproducing these results

```sh
scripts/build_decoders.sh            # patched PyMatching, Relay-BP and ldpc into the active virtualenv
scripts/pipeline.sh 10 relay         # Relay-BP data collection with 10 workers (resumable)
PYTHONPATH=. python3 scripts/analyze.py [--skip-lindley]
```

Outputs (§4–5):
- `results/exp1_relay.{json,csv}`: capped LER for every cap and point
- `results/exp2_tail_fits.{json,csv}`: all fits and tests
- `results/exp2_buffers.json`: buffer tables for $\varepsilon=10^{-6}$ to $10^{-15}$
- `results/exp2_lindley.json`: Lindley simulation
- `results/summary.json`: cap definitions
- `figures/exp1_ler_vs_p_relay.png`, `figures/exp2_ccdf_relay_*.png`

MWPM and union-find (§7):

```sh
scripts/pipeline.sh 14 mwpm uf       # data/mwpm and data/uf: d in {11,15}, p in {1e-3, 3e-3}, 1e6 shots each
PYTHONPATH=. python3 scripts/analyze_uf.py
```

Outputs: `results/uf_summary.json`, `results/uf_tail_fits.{json,csv}`, `figures/uf_ccdf_vs_mwpm.png`,
`figures/uf_components_{mwpm,uf}.png`, `figures/uf_ccdf_{mwpm,uf}_*.png`.

## 7. Union-find decoder (ldpc) vs MWPM: full operation counts

*Added 2026-10-10. Tail class only. No capping experiment was run, and no per-operation cost constants (or buffer sizes in absolute units) are estimated.*

### 7.1 Setup

- **Work model.** Every counter counts operations of bounded cost: O(1), or bounded by the fixed maximum degree of the decoding graph. Decode time is then within constant factors of the counter sum $S$. Because $c_{\min}\max_i N_i \le T \le k\,c_{\max}\max_i N_i$ holds for every decode, whatever the dependence between components, the tail class of $S$ (and of the time $T$) is that of its heaviest-tailed component, whatever the per-counter costs. So we fit the tail of $S$ and of each component group.
- **MWPM.** PyMatching (Sparse Blossom) with the extended op-counter patch (`third_party/pymatching_op_counters.patch`). There are three component groups:
  - **dual:** priority-queue pops (including stale events), radix-heap refill moves, flood events, node reschedules and their neighbour scans;
  - **primal:** matcher events, alternating-tree operations (rotations, common-ancestor walks, shattering trees into matches, path pruning), blossom-cycle and blossom-recursion region operations, and detector-node visits caused by growth-rate changes and blossom (un)wrapping;
  - **extraction:** recursive blossom shattering at the end.

  Predictions are bit-identical to stock PyMatching. The data are in `data/mwpm`. The earlier MWPM metric ($\text{pq\_pops}+\text{node\_reschedules}$, used by the removed MWPM runs) is also stored, as `S_legacy`.
- **UF.** `UnionFindDecoder` from Roffe's `ldpc` (pinned commit `d342996`) in graph (peeling) mode with unweighted growth, which adds all neighbouring bits per step. Validity is the O(1) parity + boundary test. The patch (`third_party/ldpc_uf_op_counters.patch`) replaces upstream's spanning-tree peeler, which aborts on circuit-level surface-code clusters, with the standard peeler on a BFS spanning tree; growth, merging and validity are unchanged. Corrections were checked to reproduce the syndrome exactly. There are four component groups:
  - **growth:** cluster creation, row scans, bit additions and their column visits;
  - **merge:** merges and the nodes copied (ldpc merges by copying hash sets);
  - **validity:** per-round O(1) checks and the size sort;
  - **peel:** spanning-tree construction and peeling.

  Hypergraph (matrix-inversion) UF is out of scope. The data are in `data/uf`.
- **Decoding graph and points.** Both decoders use the `decompose_errors=True` DEM of Stim's rotated memory-X circuit with $d$ rounds and uniform circuit-level depolarizing noise. Points are $d\in\{11,15\}$ and $p\in\{10^{-3},3\times10^{-3}\}$, with $10^6$ shots each and the same seeds for both decoders. Collection took about 2 minutes in total with 14 workers.

### 7.2 Results

| $d$ | $p$ | UF LER | MWPM LER | MWPM: median $S$, max/median | MWPM share of mean (dual / primal / extraction) | UF: median $S$, max/median | UF share of mean (growth / merge / validity / peel) |
|---|---|---|---|---|---|---|---|
| 11 | 1e-3 | $1.3\times10^{-5}$ | $0$ (0/10⁶) | 935, 10.9 | 0.88 / 0.10 / 0.03 | 2594, 6.2 | 0.20 / 0.05 / 0.01 / 0.74 |
| 15 | 1e-3 | $0$ (0/10⁶) | $0$ (0/10⁶) | 2644, 7.7 | 0.88 / 0.10 / 0.02 | 7310, 3.9 | 0.20 / 0.05 / 0.01 / 0.75 |
| 11 | 3e-3 | $1.0\times10^{-2}$ | $2.9\times10^{-4}$ | 3408, 12.5 | 0.88 / 0.10 / 0.02 | 8313, 4.3 | 0.19 / 0.14 / 0.01 / 0.66 |
| 15 | 3e-3 | $7.9\times10^{-3}$ | $3.7\times10^{-5}$ | 9724, 10.4 | 0.89 / 0.09 / 0.02 | 23977, 3.8 | 0.18 / 0.16 / 0.01 / 0.64 |

The two decoders' $S$ values are in different units and are not comparable with each other; only their shapes are.

Tail fits use the same pipeline as §5 (threshold at the top $10^{-3}$, about 1000 exceedances). Each entry is the GPD shape $\xi$ with its 95% bootstrap CI, then in brackets $\xi$ at a 10× higher threshold.

| $d$ | $p$ | MWPM $S$ | MWPM dual | MWPM primal | MWPM extraction |
|---|---|---|---|---|---|
| 11 | 1e-3 | 0.23 (0.16, 0.28) [0.26] | 0.20 (0.15, 0.26) [0.21] | 0.21 (0.17, 0.26) [0.15] | −0.07 (−0.10, −0.03) |
| 15 | 1e-3 | 0.17 (0.12, 0.22) [0.20] | 0.17 (0.12, 0.22) [0.26] | 0.24 (0.17, 0.28) [0.26] | −0.04 (−0.10, 0.01) |
| 11 | 3e-3 | 0.23 (0.18, 0.28) [0.18] | 0.20 (0.16, 0.26) [0.13] | **0.37** (0.33, 0.44) [0.30] | −0.02 (−0.06, 0.02) |
| 15 | 3e-3 | 0.31 (0.28, 0.36) [0.18] | 0.27 (0.22, 0.33) [0.12] | **0.45** (0.39, 0.53) [0.25] | −0.02 (−0.06, 0.01) |

| $d$ | $p$ | UF $S$ | UF growth | UF merge | UF validity | UF peel |
|---|---|---|---|---|---|---|
| 11 | 1e-3 | 0.15 (0.10, 0.19) [0.09] | 0.13 (0.07, 0.18) [0.06] | 0.08 (0.03, 0.14) | −0.13 | 0.13 (0.09, 0.19) [0.04] |
| 15 | 1e-3 | 0.14 (0.09, 0.19) [0.26] | 0.17 (0.11, 0.23) [0.27] | 0.06 (0.02, 0.10) | 0.00 | 0.10 (0.06, 0.14) [0.27] |
| 11 | 3e-3 | −0.01 (−0.05, 0.03) [−0.20] | 0.01 (−0.04, 0.06) [−0.09] | 0.02 (−0.02, 0.06) | 0.04 | −0.03 (−0.07, 0.01) [−0.18] |
| 15 | 3e-3 | 0.04 (−0.00, 0.08) [−0.05] | 0.10 (0.06, 0.14) [0.09] | 0.05 (0.00, 0.09) | −0.17 | 0.04 (−0.00, 0.08) [−0.09] |

The best AIC model for $S$ is Pareto at all four MWPM points and at the two UF points with $p=10^{-3}$. At the UF points with $p=3\times10^{-3}$ it is exponential ($d=11$) and lognormal ($d=15$). All fits, Vuong tests and per-counter means and maxima are in `results/uf_tail_fits.csv` and `results/uf_summary.json`.

Figures:
- `figures/uf_ccdf_vs_mwpm.png`: survival functions of $S$, normalized by the median;
- `figures/uf_components_{mwpm,uf}.png`: each component group against the total;
- `figures/uf_ccdf_{mwpm,uf}_*.png`: empirical log-log survival function of $S$ for each point (no fitted curves).

### 7.3 Findings

- **MWPM's heaviest tail is on the primal side.** Primal work (alternating trees, blossoms and the region-wide updates they trigger) is only about 10% of the mean but has the heaviest tail. At $p=3\times10^{-3}$, its $\xi=0.37$–$0.45$ (Pareto $\alpha=1/\xi\approx2.2$–$2.7$) against 0.20–0.27 for the dual. Its tail also gets heavier with $p$ and $d$. The ordering persists at the higher threshold (0.25–0.30 against 0.12–0.13). Within the observed range ($\ge10^{-6}$), $S$ is dominated by the dual, so the fit to $S$ itself gives a smaller $\xi$ (0.17–0.31). Asymptotically, by the bound above, the primal component sets the tail class of MWPM work.
- **The earlier MWPM metric (`S_legacy`) understates MWPM's tail.** It counts only PQ pops and reschedules. It sees primal events only through the reschedules they trigger, and misses the alternating-tree bookkeeping, blossom-cycle work and node visits of shrinking regions. It is mostly a dual metric. For capping and buffer conclusions at large $d$ and $p$, the primal tail should be used instead.
- **UF's tail is lighter than MWPM's at every point, and lightens with $p$.**
  - The worst of $10^6$ UF decodes is 3.8–6.2× the median, against 7.7–12.5× for MWPM.
  - At $p=10^{-3}$, UF's $\xi\approx0.14$–$0.15$ (CI above 0, so a mildly heavy, Pareto-type tail with $\alpha\approx7$).
  - At $p=3\times10^{-3}$, $\xi\approx0$ (CI includes 0; Gumbel domain).
  - The opposite $p$ trends make sense. UF's work is a sum over many fairly independent clusters, which smooths out as the number of defects grows, while MWPM's primal work concentrates in large alternating trees and blossoms, which become more common as $p$ grows.
  - Among UF's components, growth (cluster-boundary scans) usually has the largest $\xi$. Peeling dominates the mean (64–75%). The higher-threshold fits move UF's $\xi$ by up to about 0.1–0.2 (for example, $S$ at $d=15$, $p=10^{-3}$: 0.14 → 0.26), so the ordering of UF's components is not well resolved. The robust statement is $\xi\lesssim0.25$ at $p=10^{-3}$ and $\xi\approx0$ at $p=3\times10^{-3}$.
- **This supersedes the earlier cluster-operation result.** The earlier version of this section found a bounded-type tail ($\xi\approx-0.1$). That metric charged one unit per cluster-level event, regardless of cluster size. It was a light-tailed artifact of following the detection-event count.
- **Accuracy caveat.** This UF variant is far less accurate than MWPM under circuit-level noise: about 36× ($d=11$) and 210× ($d=15$) higher LER at $p=3\times10^{-3}$. Contributing causes are unweighted growth by full edges and peeling corrections that are not minimum-weight. $p=3\times10^{-3}$ is close to this decoder's threshold: its LER only falls from $1.0\times10^{-2}$ to $7.9\times10^{-3}$ between $d=11$ and $d=15$.
- **Statistics caveat.** With $10^6$ shots the fits describe the distribution down to about $10^{-6}$. Statements about tail classes beyond that level assume the fitted class continues.
