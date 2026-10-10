# synbuf: Results Report (Experiments 1 and 2)

*Status as of 2026-10-04. Describes the runs actually performed and the results currently on disk in `data/`, `results/` and `figures/`. The plan these runs implement is in [EXPERIMENTS.md](EXPERIMENTS.md); deviations from it are listed in [§2](#2-what-was-and-was-not-run).*

## 1. Summary of findings

1. **Capping the decoder blows up the logical error rate; it does not, in general, produce a flat error floor.** With a fixed work budget $c$, the capped LER $P_L(c)$ tracks the uncapped LER at low $p$ and then rises far faster than it as $p$ grows, eventually saturating at 1. The excess is set by $P(S>c)$, which is strongly $p$- and $d$-dependent, so the cap acts as an additional, rapidly growing failure mechanism rather than a constant offset.
2. **The blow-up gets worse with distance.** Decode work grows with $d$ (MWPM median work at $p=10^{-3}$: 14, 62 and 167 operations at $d=7, 11, 15$), so a fixed budget binds harder on larger codes. At $p=3\times10^{-3}$ with $c=1000$, the MWPM LER is $1.7\times10^{-3}$ at $d=7$, $4.5\times10^{-4}$ at $d=11$, and $1.4\times10^{-2}$ at $d=15$: going from $d=11$ to $d=15$ makes the capped decoder about 30× *worse*, while the uncapped LER improves about 8×. Capping is therefore incompatible with scaling distance toward $\varepsilon\approx10^{-12}$, which **motivates buffering syndromes instead of capping**.
3. **MWPM work has a mildly heavy tail.** At $d=7$ the tail is close to the Gumbel domain (lognormal-like, GPD shape $\xi\approx0.04$–$0.12$). At $d=11$ and $15$ it is Pareto-type with $\xi\approx0.14$–$0.28$ (Pareto $\alpha\approx5.7$–$9.7$). The fitted tails imply modest buffers: tens to a few hundred median decodes for a sliding window at $\varepsilon=10^{-12}$, and under 15 median decodes with 2 spare servers.
4. **Relay-BP work has a much heavier, power-law tail that grows heavier with $p$.** $\xi$ ranges from about 0.05–0.4 at $p=5\times10^{-4}$, through 0.34–0.94 at $10^{-3}$, to above 1 at $3\times10^{-3}$ for the three codes with $d\ge10$ ($\alpha\approx0.74$–$1.01$). A fit at a 10× higher threshold gives $\xi\approx0.84$–$1.15$ there, so the tail sits at or near the infinite-mean boundary rather than clearly beyond it. With $\alpha\le1$ the leading-order buffer formulas diverge, and a D/G/1 Lindley simulation of [[144,12,12]] at $p=3\times10^{-3}$ spends about a quarter of its time draining single huge decodes even at load $\rho=0.5$. Spare servers (parallel windows) shrink the Pareto-tail buffer by orders of magnitude where the formulas still apply.

5. **Union-find (ldpc) work, counted in cluster operations, is light-tailed** ([§7](#7-union-find-decoder-ldpc-work-distribution)). At $d\in\{11,15\}$ and $p\in\{10^{-3},3\times10^{-3}\}$, the GPD shape is $\xi\approx-0.1$ (bounded-type tail), the maximum is at most 3.3× the median, and the $\varepsilon=10^{-12}$ buffer is 2–4 median decodes. Counted per node touched, ldpc's merge-by-copy makes the tail Gumbel-like (max up to 13× the median). This UF variant is much less accurate than MWPM: about 34× and 200× higher LER at $p=3\times10^{-3}$ for $d=11$ and $d=15$.

## 2. What was and was not run

### Run

| Item | What was done |
|---|---|
| MWPM data | Rotated surface code memory, $d\in\{3,7,11,15\}$, $p\in\{5\times10^{-4},10^{-3},3\times10^{-3}\}$, $10^6$ shots per point (12 points) |
| Relay-BP data | BB codes [[72,12,6]], [[90,8,10]], [[108,8,10]], [[144,12,12]], same $p$ grid, $10^6$ shots per point (12 points) |
| Instrumented decoders | PyMatching and Relay-BP built from pinned commits with the patches in `third_party/` (`scripts/build_decoders.sh`) |
| Relay-BP parameters | Per-code tuned values from `results/relay_params.json` (tuned earlier; see `results/relay_tuning_notes.md`) |
| Experiment 1 | Offline $P_L(c)$ for an absolute cap ladder, the deployed Relay-BP limit, and calibrated caps $c_3,c_4,c_5$ |
| Experiment 2 | Peaks-over-threshold tail fits (GPD, exponential, lognormal, Pareto), Vuong tests, AIC ranking and censored likelihood for every MWPM point with $d\ge5$ and every Relay-BP point; diagnostics and buffer tables for 9 selected points; extrapolation plots for every point with $d\ge7$; D/G/1 Lindley simulation at 2 points |

### Not run, or run differently from EXPERIMENTS.md

| Planned in EXPERIMENTS.md | Status |
|---|---|
| Odd $d$ from 3 to about 15 | Only $d=3,7,11,15$; $d=5,9,13$ were skipped. $d=3$ is excluded from the tail fits by design. |
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
| Lindley check | D/G/1 only, at 2 points; numbers in `results/exp2_lindley.json`, figure removed |

## 3. Common setup

**Circuits.** Generated with Stim (`synbuf/circuits.py`).
- Surface code: Stim's rotated memory circuit with $d$ rounds, with `after_clifford_depolarization`, `before_round_data_depolarization` and `before_measure_flip_probability` all set to $p$.
- BB codes: the Bravyi et al. (Nature 627, 778, 2024) syndrome cycle with uniform circuit-level noise, $d$ rounds.

**Window.** One $d$-round memory block is decoded as one window. All LERs are per window, and for BB codes per window of the whole block (any logical error among the $k$ logical qubits counts as one failure).

**Work metric $S$.**
- MWPM: $S$ = priority-queue pops + detector-node reschedules in Sparse Blossom, from the op-counter patch (`decode_batch_with_op_counts`). Flooder and matcher event counts are recorded too but not used.
- Relay-BP: $S$ = total BP iterations up to and including the first converged leg (time to first solution, `stop_nconv=1`). Each decode runs 80 pre-iterations, then up to 10,000 legs of up to 60 iterations each, so the leg limit is $80+10{,}000\times60=600{,}080$ iterations. A shot that reaches it without converging is **censored** and counted as a logical failure.

**Sampling.** Shots are split into chunks with per-chunk seeds and decoded by 10 worker processes, so results do not depend on the worker count. Data: `data/mwpm/surface_d{d}_p{p}.npz` and `data/relay/bb{code}_p{p}/`.

## 4. Experiment 1: decoder caps and the LER

### 4.1 Method

Every shot is decoded to completion and records $(S, \mathbb 1[F_\infty])$. For any cap $c$ the capped LER is computed offline, with an aborted decode counted as a failure:

$$P_L(c) = P(F_\infty,\ S\le c) + P(S>c).$$

### 4.2 Caps

- **Absolute ladder.** Caps on a 1-2-5 grid from the smallest median $S$ to the largest observed $S$, kept below the Relay-BP leg limit: 2 to 2,000 operations for MWPM, 2 to 500,000 iterations for Relay-BP. Each cap is held fixed while $d$ and $p$ vary.
- **Deployed Relay-BP limit.** $80 + 300\times60 = 18{,}080$ iterations (300 legs, the limit used during tuning).
- **Calibrated caps.** $c_j$ with $P(S>c_j)=10^{-j}$ at a reference point. For MWPM the reference is $(d=11, p=10^{-3})$, giving $c_3=209$, $c_4=295$, $c_5=442$. For Relay-BP the reference is ([[144,12,12]], $p=3\times10^{-3}$), giving $c_3=2{,}806$ and $c_4=66{,}456$; $c_5$ cannot be calibrated because the censored fraction there ($2.6\times10^{-5}$) exceeds $10^{-5}$. These are listed in `results/exp1_*.{json,csv}` but not plotted.

### 4.3 Results

Figures: [`figures/exp1_ler_vs_p_mwpm.png`](figures/exp1_ler_vs_p_mwpm.png), [`figures/exp1_ler_vs_p_relay.png`](figures/exp1_ler_vs_p_relay.png). Each panel is one code. It shows the uncapped LER (dashed black) and up to five ladder caps that bind somewhere on the $p$ grid but cap fewer than 5% of shots at the lowest $p$; for Relay-BP it also shows the deployed limit (purple).

**MWPM, $P_L(c)$ per window** (a value of 0 means no failures in $10^6$ shots):

| $d$ | $p$ | median $S$ | uncapped | $c=200$ | $c=500$ | $c=1000$ | $c=2000$ |
|---|---|---|---|---|---|---|---|
| 7 | 5e-4 | 6 | 1.0e-6 | 3.0e-6 | 1.0e-6 | 1.0e-6 | 1.0e-6 |
| 7 | 1e-3 | 14 | 2.0e-5 | 3.8e-5 | 2.0e-5 | 2.0e-5 | 2.0e-5 |
| 7 | 3e-3 | 49 | 1.71e-3 | 6.18e-3 | 1.72e-3 | 1.71e-3 | 1.71e-3 |
| 11 | 5e-4 | 30 | 0 | 3.1e-5 | 0 | 0 | 0 |
| 11 | 1e-3 | 62 | 0 | 1.35e-3 | 6.0e-6 | 0 | 0 |
| 11 | 3e-3 | 225 | 2.9e-4 | 0.648 | 9.23e-3 | 4.52e-4 | 2.91e-4 |
| 15 | 5e-4 | 79 | 0 | 1.08e-3 | 1.0e-6 | 0 | 0 |
| 15 | 1e-3 | 167 | 0 | 0.230 | 7.1e-5 | 1.0e-6 | 0 |
| 15 | 3e-3 | 612 | 3.7e-5 | 1.00 | 0.848 | 1.43e-2 | 1.54e-4 |

**Relay-BP, $P_L(c)$ per window** (deployed limit = 18,080 iterations):

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

- **Blow-up, not a floor.** For a fixed cap, $P_L(c)$ converges to the uncapped LER as $p$ falls and rises much more steeply than it as $p$ grows. For MWPM at $d=15$, $c=1000$: below $10^{-6}$ at $p=5\times10^{-4}$, $10^{-6}$ at $10^{-3}$, and $1.4\times10^{-2}$ at $3\times10^{-3}$ (390× the uncapped LER). A floor would require $P(S>c)$ to be insensitive to $p$, but it is not: the cap opens a new failure channel whose rate grows quickly with $p$ and $d$.
- **Larger distance makes capping worse.** At fixed $p$ and $c$, $P_L(c)$ is non-monotone in $d$ or increases with it. At $p=10^{-3}$, $c=500$: $2.0\times10^{-5}$, $6\times10^{-6}$ and $7.1\times10^{-5}$ at $d=7, 11, 15$, while the uncapped LER at $d=11$ and $15$ is below $10^{-6}$. The cap that is harmless at $d=7$ dominates at $d=15$.
- **Relay-BP needs caps far above its median.** Because of its heavy tail, caps of $10^2$–$10^4\times$ the median still add measurable error at $p\ge10^{-3}$ for the $d\ge10$ codes. The deployed 300-leg limit raises the LER by 1.3–1.5× at $p=3\times10^{-3}$ and by about 3–6× at $p=10^{-3}$ for [[90,8,10]] and [[108,8,10]], though the low-$p$ counts are small (4 uncapped vs 24 capped failures for [[90,8,10]]). [[72,12,6]] is insensitive to every cap of 2,000 iterations or more because its LER is high.
- **Statistics.** At $p\le10^{-3}$ for $d\ge11$ (MWPM) and the $d\ge10$ BB codes, the uncapped LER is at or below the $10^{-6}$ resolution of $10^6$ shots, so ratios at those points are lower bounds or rest on a handful of events. Error bars in the figures are 1σ Wilson intervals.

## 5. Experiment 2: tail class of the work distribution

### 5.1 Method

Implemented in `synbuf/tails.py` and `scripts/analyze.py`, on the uncapped Experiment 1 data.
- **Threshold.** $u$ is the empirical $(1-q)$ quantile of $S$, with $q=10^{-3}$ for MWPM and $q=10^{-2}$ for Relay-BP (about 2,000 and 10,000 exceedances respectively). A second fit at $q/10$ checks threshold stability.
- **Fits.** Left-truncated MLE fits above $u$ of exponential, lognormal and Pareto tails, plus a generalized Pareto (GPD) fit with a bootstrap 95% interval on its shape $\xi$. $\xi>0$ indicates a Pareto-type tail; $\xi\approx0$ the Gumbel domain (exponential or lognormal).
- **Model selection.** AIC ranking and pairwise Vuong likelihood-ratio tests.
- **Censoring.** Relay-BP shots at the leg limit enter the likelihood as right-censored observations ($S>600{,}080$).
- **Buffers.** The fitted tail is spliced onto the empirical bulk below $u$. For sliding windows (D/G/1), solve $\frac{\rho}{1-\rho}\bar F_I(x)=\varepsilon$ (Veraverbeke), where $\bar F_I$ is the integrated tail. For parallel windows with $m$ spare servers, solve $\bar F_I(x)^m=\varepsilon$ (Foss–Korshunov, leading order). $\rho=0.5$.
- **Lindley check.** A D/G/1 queue simulated for $2\times10^8$ windows by resampling the empirical $S$, at $\rho=0.5$ and $0.8$.

### 5.2 Figures

- Extrapolation of each fitted tail to $P(S>x)=10^{-15}$, for every point with $d\ge7$ (18 plots): `figures/exp2_extrap_{mwpm,relay}_*.png`.
- Three-panel diagnostics ($\log\bar F$ vs $x$, $\log\bar F$ vs $\log x$ with fits, normal Q–Q of $\log S$) for 9 selected points: `figures/exp2_diag_*.png`.

`figures/exp2_param_trends.png`, `exp2_buffer_vs_eps.png` and `exp2_lindley.png` are still on disk from an earlier run, but the analysis no longer generates them.

### 5.3 Tail fits

From `results/exp2_tail_fits.csv`. "Censored" is the number of shots at the leg limit. "Vuong LN|P $p$" is the p-value of the lognormal-vs-Pareto test.

**MWPM:**

| $d$ | $p$ | $u$ | Best (AIC) | $\xi$ [95% CI] | Pareto $\alpha$ | Vuong LN\|P $p$ |
|---|---|---|---|---|---|---|
| 7 | 5e-4 | 53 | lognormal | 0.12 [0.07, 0.15] | 4.6 | 0.05 |
| 7 | 1e-3 | 87 | lognormal | 0.05 [0.01, 0.10] | 4.9 | <0.01 |
| 7 | 3e-3 | 240 | lognormal | 0.04 [−0.01, 0.08] | 5.4 | <0.01 |
| 11 | 5e-4 | 106 | Pareto | 0.14 [0.09, 0.18] | 6.6 | 0.98 |
| 11 | 1e-3 | 190 | Pareto | 0.18 [0.12, 0.22] | 6.8 | 0.69 |
| 11 | 3e-3 | 640 | Pareto | 0.21 [0.17, 0.26] | 5.7 | 0.65 |
| 15 | 5e-4 | 188 | Pareto | 0.14 [0.09, 0.18] | 9.4 | 0.57 |
| 15 | 1e-3 | 351 | Pareto | 0.14 [0.07, 0.19] | 9.7 | 0.76 |
| 15 | 3e-3 | 1261 | Pareto | 0.28 [0.22, 0.32] | 6.5 | 0.76 |

**Relay-BP:**

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
| MWPM $d=15$, $p=5\times10^{-4}$ | Pareto | 26 | 5.0 | 2.9 | 48 |
| MWPM $d=15$, $p=10^{-3}$ | Pareto | 21 | 4.3 | 2.5 | 40 |
| MWPM $d=11$, $p=3\times10^{-3}$ | Pareto | 248 | 12.9 | 4.8 | 547 |
| MWPM $d=15$, $p=3\times10^{-3}$ | Pareto | 87 | 6.9 | 3.0 | 2,350 |
| Relay [[144,12,12]], $p=10^{-3}$ | Pareto | $9.5\times10^4$ | 247 | 34 | $4\times10^9$ |
| Relay [[72,12,6]], $p=3\times10^{-3}$ | lognormal | $5.1\times10^7$ | $1.2\times10^5$ | $8.0\times10^3$ | $5.3\times10^{13}$ |
| Relay [[90,8,10]], [[108,8,10]], [[144,12,12]], $p=3\times10^{-3}$ | Pareto | diverges | diverges | diverges | diverges |

### 5.5 Lindley check (D/G/1)

From `results/exp2_lindley.json`. $P(W>x)$ is the fraction of windows whose backlog exceeds $x$, with $x$ in multiples of the median $S$.

| Point | $\rho$ | $x=5$ | $x=50$ | $x=500$ | $x=5000$ |
|---|---|---|---|---|---|
| MWPM $d=15$, $p=3\times10^{-3}$ | 0.5 | $1.3\times10^{-6}$ | 0 | 0 | 0 |
| MWPM $d=15$, $p=3\times10^{-3}$ | 0.8 | $1.4\times10^{-5}$ | 0 | 0 | 0 |
| Relay [[144,12,12]], $p=3\times10^{-3}$ | 0.5 | 0.40 | 0.38 | 0.34 | 0.25 |
| Relay [[144,12,12]], $p=3\times10^{-3}$ | 0.8 | 0.74 | 0.71 | 0.68 | 0.59 |

For Relay-BP, each decode that reaches the leg limit (about 1 in 40,000 windows) takes about 8,000 windows to drain, so the queue spends about a quarter of its time recovering from one such decode. This is the single-big-jump regime in its extreme form. The Veraverbeke asymptotic exceeds 1 at these thresholds and is not valid here.

### 5.6 Findings

- **MWPM.**
  - The bulk-lognormal hypothesis holds at $d=7$, where the tail is near the Gumbel domain and lognormal beats Pareto in the Vuong test.
  - At $d\ge11$ the far tail is Pareto-type ($\xi\approx0.14$–$0.28$, CIs excluding 0). Exponential is rejected throughout. Lognormal and Pareto cannot be distinguished: the lognormal fit degenerates to a power law (very large $\sigma$, very negative $\mu$) and loses only on AIC's parameter penalty.
  - Because $\alpha\approx6$–$10$, buffers stay modest: 21–248 median decodes for a sliding window at $10^{-12}$, and 2.5–13 with 2–3 spare servers. The GPD fits are more pessimistic (up to 2,350 at $d=15$, $p=3\times10^{-3}$), so the far-tail extrapolation is the main uncertainty.
- **Relay-BP.**
  - The tail is heavier than MWPM's at every $p\ge10^{-3}$ and grows heavier with $p$.
  - At $p=3\times10^{-3}$ the three $d\ge10$ codes have $\alpha\le1$ over the observed range. The single-server buffer has no finite leading-order size, and the simulated queue is dominated by individual huge decodes.
  - Where $\alpha>1$ ([[144,12,12]] at $p=10^{-3}$), spare servers cut the buffer from about $10^5$ to $34$ median decodes ($m=3$). This matches the prediction that parallel windows help most for a fat tail.
- **Threshold caveat.** For Relay-BP the threshold $u$ sits only a few times above the median at $p\le10^{-3}$ (11–24 iterations), so those fits describe the top 1% of the distribution rather than only the far tail. The higher-threshold fits in `exp2_tail_fits.csv` (`xi_high_u`, `best_aic_high_u`) also give $\xi>0$ at every point but move $\xi$ by up to about 1 at some points (for example [[108,8,10]] at $p=5\times10^{-4}$: 0.42 vs 1.46). At $p=3\times10^{-3}$ the higher threshold gives $\xi=0.86$, $0.84$ and $1.15$ for [[90,8,10]], [[108,8,10]] and [[144,12,12]], and prefers lognormal for the first two. So $\xi>1$ at that $p$ is threshold-sensitive; the robust statement is that the tail is at or near the infinite-mean boundary.
- **Censoring caveat.** Relay-BP censoring is at most $1.3\times10^{-4}$ of shots. It affects only the tail beyond 600,080 iterations, which is irrelevant for any practical cap but means buffer estimates that extrapolate past that point are effectively bounded by it.

## 6. Reproducing these results

```sh
scripts/build_decoders.sh            # patched PyMatching + Relay-BP into the active virtualenv
scripts/pipeline.sh 10 mwpm relay    # data collection with 10 workers (resumable)
PYTHONPATH=. python3 scripts/analyze.py [--skip-lindley]
```

Outputs:
- `results/exp1_{mwpm,relay}.{json,csv}`: capped LER for every cap and point
- `results/exp2_tail_fits.{json,csv}`: all fits and tests
- `results/exp2_buffers.json`: buffer tables for $\varepsilon=10^{-6}$ to $10^{-15}$
- `results/exp2_lindley.json`: Lindley simulation
- `results/summary.json`: cap definitions
- `figures/`: the plots referenced above

Union-find (§7):

```sh
scripts/pipeline.sh 14 uf            # data/uf/, d in {11,15}, p in {1e-3, 3e-3}, 1e6 shots each
PYTHONPATH=. python3 scripts/analyze_uf.py
```

Outputs: `results/uf_summary.json`, `results/uf_tail_fits.{json,csv}`, `results/uf_buffers.json`,
`figures/uf_ccdf_vs_mwpm.png`, `figures/uf_diag_*.png`, `figures/uf_extrap_*.png`.

## 7. Union-find decoder (ldpc): work distribution

*Added 2026-10-10. Tail characterization only; no capping experiment was run for UF.*

### 7.1 Setup

- **Decoder.** `UnionFindDecoder` from Roffe's `ldpc` package (pinned commit `d342996`, see `third_party/README.md`), instrumented with per-decode cluster-operation counters (`third_party/ldpc_uf_op_counters.patch`). The algorithm is unchanged.
- **Mode.** Matrix-inversion mode with unweighted growth (all neighbouring bits per step). Each grown cluster is checked and solved by GF(2) elimination on its sub-matrix. The peeling mode was not usable: on the circuit-level graph it aborts (`Iterator index out of bounds`) on clusters that touch several boundary edges. Weighted growth (`llrs`, `bits_per_step`=1, 4, 16) gave a *higher* LER than unweighted growth in a pilot run at $d=5$, so it was not used.
- **Decoding graph.** Same as MWPM: the `decompose_errors=True` DEM of Stim's rotated memory-X circuit, $d$ rounds, uniform circuit-level depolarizing noise. Each graphlike component is a column of $H$ (weight 1 for a boundary edge).
- **Points.** $d\in\{11,15\}$, $p\in\{10^{-3},3\times10^{-3}\}$, $10^6$ shots each (`synbuf/run_uf.py`). Total collection time was about 3 minutes with 14 workers.
- **Work metric.** $S$ = number of cluster operations = singleton clusters created (one per syndrome bit) + growth steps of invalid clusters (each followed by one validity solve) + cluster merges. As a cross-check the analysis also uses a node-level cost $W$ = bits added to clusters + nodes copied by merges. ldpc merges clusters by copying hash sets rather than through a disjoint-set forest, so $W$ is closer to the implementation's actual memory traffic.

### 7.2 Results

| $d$ | $p$ | UF LER | MWPM LER | mean $S$ | median $S$ | max $S$ | max/median: UF $S$ | max/median: MWPM | max/median: UF $W$ |
|---|---|---|---|---|---|---|---|---|---|
| 11 | 1e-3 | $1.9\times10^{-5}$ | $0$ (0/10⁶) | 47.3 | 46 | 153 | 3.3 | 10.1 | 13.0 |
| 15 | 1e-3 | $0$ (0/10⁶) | $0$ (0/10⁶) | 124.8 | 124 | 286 | 2.3 | 7.3 | 6.5 |
| 11 | 3e-3 | $1.0\times10^{-2}$ | $2.9\times10^{-4}$ | 146.7 | 145 | 312 | 2.2 | 10.4 | 11.3 |
| 15 | 3e-3 | $7.4\times10^{-3}$ | $3.7\times10^{-5}$ | 388.6 | 387 | 660 | 1.7 | 8.0 | 8.4 |

Tail fits use the same pipeline and threshold as MWPM (top $10^{-3}$, about 1000 exceedances):

| $d$ | $p$ | GPD $\xi$ (95% CI) | $\xi$ at 10× higher threshold | best AIC | lognormal $\sigma$ | node cost $W$: $\xi$ (95% CI) |
|---|---|---|---|---|---|---|
| 11 | 1e-3 | −0.105 (−0.143, −0.076) | −0.095 | lognormal | 0.120 | 0.031 (−0.043, 0.093) |
| 15 | 1e-3 | −0.083 (−0.120, −0.053) | −0.052 | lognormal | 0.110 | 0.042 (−0.023, 0.096) |
| 11 | 3e-3 | −0.106 (−0.146, −0.073) | −0.124 | lognormal | 0.097 | 0.047 (−0.021, 0.102) |
| 15 | 3e-3 | −0.104 (−0.149, −0.067) | −0.048 | lognormal | 0.068 | 0.042 (−0.028, 0.109) |

Buffer at $\varepsilon=10^{-12}$ in median decodes ($\rho=0.5$, lognormal tail spliced onto the empirical bulk, `results/uf_buffers.json`): 4.2 / 2.9 / 2.7 / 1.9 for D/G/1 at ($d$,$p$) = (11,1e-3) / (15,1e-3) / (11,3e-3) / (15,3e-3). With two spare servers these fall to 3.0 / 2.1 / 2.0 / 1.6. The exponential fit gives slightly larger values (5.5 / 3.5 / 3.2 / 2.3 for D/G/1).

Figures: `figures/uf_ccdf_vs_mwpm.png` (UF and MWPM survival functions, normalized by their medians), `figures/uf_diag_*.png` (the three straight-line diagnostics) and `figures/uf_extrap_*.png` (tail fits extrapolated to $10^{-15}$).

### 7.3 Findings

- **The cluster-operation count has a light, bounded-type tail.** $\xi<0$ with the whole 95% CI below 0 at all four points (Weibull domain). The worst of $10^6$ decodes is only 1.7–3.3× the median, against 7–10× for MWPM's operation count at the same points. The tail also gets *relatively* lighter as $d$ and $p$ grow. $S$ has a rank correlation of 0.99 with the number of syndrome bits, so its distribution essentially follows the detection-event count, which is a sum of many independent Bernoulli variables. The resulting buffer is only a few median decodes at $\varepsilon=10^{-12}$.
- **At node granularity the tail is heavier.** $W$ (bits added + nodes copied by merges) reaches 6.5–13× its median, with $\xi\approx0.03$–$0.05$ and a CI that includes 0 (Gumbel domain, exponential or lognormal). That is comparable to MWPM at $d=7$. The extra spread comes from merges of large clusters, which copy their node sets in this implementation. A disjoint-set-forest UF with union by size would bound that cost. So the light tail is a property of the cluster-level operation count, and of a UF implementation whose merges are near-constant time, not of ldpc's merge-by-copy.
- **Accuracy caveat.** This UF variant is far less accurate than MWPM under circuit-level noise: 34× (d=11) and 200× (d=15) higher LER at $p=3\times10^{-3}$. Part of this is unweighted growth and part is that inversion mode returns an arbitrary rather than a minimum-weight correction within each cluster. $p=3\times10^{-3}$ is close to this decoder's threshold: the LER only falls from $1.0\times10^{-2}$ to $7.4\times10^{-3}$ between $d=11$ and $d=15$.
- **Statistics caveat.** With $10^6$ shots the fits describe the distribution down to about $10^{-6}$. The extrapolations to $10^{-12}$ in the buffer table assume the fitted class continues beyond that.
