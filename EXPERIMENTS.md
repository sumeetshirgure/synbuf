# synbuf — Experiment Plan

*Syndrome buffer sizing for real-time QEC decoding. Target: a short letter (IEEE Computer Architecture Letters as first choice).*

## 1. The problem in one paragraph

A quantum memory emits a deterministic stream of syndrome bits. The decoder's work per window is random, so a classical buffer has to absorb slow decodes. Two options exist:

- **Cap the decoder.** Abort decodes that run too long. Aborted windows become logical failures, which can create an LER floor.
- **Don't cap.** Then the buffer must be large enough that its overflow probability stays below the target logical error rate $\varepsilon \approx 10^{-12}$.

How large depends almost entirely on the **tail class** of the decoder's work distribution. synbuf measures that tail class for MWPM and Relay-BP, shows when capping causes an error floor, and turns the tail class into buffer-size scaling laws for sliding-window (D/G/1) and parallel-window (D/G/k) decoding.

## 2. Notation

| Symbol | Meaning |
|---|---|
| $S$ | Algorithmic work for one window decode (random variable). Measured as an operation count, not wall-clock time. |
| $c$ | Cap: work limit after which a decode is aborted |
| $a$ | Commit period of one window, in the same work units (work budget per committed window) |
| $\rho = E[S]/a$ | Load (single server); $\rho < k$ required for $k$ servers |
| $k$ | Number of decoder servers |
| $m \approx \lceil k-\rho \rceil$ | Number of spare servers: how many can be stuck before the system overloads |
| $F_\infty$ | Event: logical failure of the window when decoded with no cap |
| $F_c$ | Event: logical failure of the window with cap $c$ (an aborted decode counts as a failure) |
| $p_L^{(0)} = P(F_\infty)$ | Uncapped logical error rate per window |
| $P_L(c) = P(F_c)$ | Capped logical error rate per window |
| $\varepsilon$ | Target total failure probability per window (about $10^{-12}$) |
| $d,\ p$ | Code distance; physical (circuit-level) error rate |

**Unit of accounting:** always per window (or consistently per $d$ rounds). Never mix per-shot, per-round and per-window rates.

## 3. Scope and modeling assumptions

**Decoder and code pairs**
- MWPM (PyMatching / Sparse Blossom) on the rotated surface code.
- Relay-BP on bivariate bicycle (BB) codes, with its practical iteration and leg limits relaxed so the true tail is visible.

**Explicitly out of scope:** parallel-window decoding for BP (a separate decoder contribution), fallback decoders, and correlations between windows.

**Assumptions to state in the paper**
1. **Independence.** Under i.i.d. circuit noise, disjoint windows have essentially disjoint fault sets, so window work is modeled as i.i.d. Overlap correlations would change constants only.
2. **Sliding window = D/G/1.** Arrivals are deterministic (a window is ready when its last round lands), and there is one server because window $i+1$ needs window $i$'s commits.
3. **Parallel windows ≈ D/G/k.** Capacity is set by $\lambda E[S_A + S_B] < k$, where $S_A$ and $S_B$ are the A- and B-window work.
4. **Out-of-order memory release** (block pool). This is load-bearing. With an in-order ring buffer, one stuck window pins all later memory, the tail stays single-jump, and extra servers do not shrink the buffer.
5. **Stability** ($\rho < 1$ or $\rho < k$). The throughput argument for parallel windows is kept separate from the tail argument.
6. **Decoder switching** (Toshio et al., arXiv:2510.25222). Under double-window decoding the weak decoder's queue is D/G/1 regardless of the strong decoder. The weak decoder's service time must include computing its soft output; the complementary gap needs a second matching.

## 4. Experiment 1: Cap vs LER floor

**Goal.** Show that a fixed decoder work budget $c$ creates an LER floor, and that the floor gets *worse* with distance. That makes capping incompatible with scaling to $\varepsilon \approx 10^{-12}$ and motivates the buffer.

### Setup
- **Rotated surface code + MWPM:** odd $d$ from 3 up to about 15. Circuit-level noise (uniform depolarizing or SI1000). $p \in \{5\times10^{-4},\ 10^{-3},\ 3\times10^{-3}\}$.
- **BB codes + Relay-BP:** same $p$ grid as MWPM; [[72,12,6]], [[90,8,10]], [[108,8,10]], [[144,12,12]] (optionally coprime BB codes for more points). Use the Bravyi et al. syndrome circuit in Stim and validate LERs against published values first.
  - Caveat: these codes are not a smooth family; $n$ and $k$ change irregularly with $d$. Use the two $d=10$ codes as a sanity check on structure versus distance.
  - Retune Relay-BP parameters per code. They were tuned for the gross code.
  - Normalize LER per window of the whole block, since $k$ varies.
- **Windows:** a $d$-round memory block decoded as one window is enough for the letter.

### Work metric (not wall-clock time)
- **MWPM:** instrument Sparse Blossom to count elementary operations (e.g. flooder events / priority-queue pops). This requires a small fork of PyMatching. Wall-clock time on a workstation is dominated by OS jitter in the tail. If wall-clock is used as a stopgap, pin the core and treat it as an upper bound.
- **Relay-BP:** total BP iterations summed over legs. Set the leg limit effectively infinite. Report **time to first solution** (stopping after $S=1$ solution) separately from time to collect $S>1$ solutions, because the latter is a design choice, not tail.

### Procedure
For every shot, record the pair $(S,\ \mathbb{1}[F_\infty])$, decoding to completion. Then, for any cap $c$, compute offline (no re-decoding needed):

$$P_L(c) = P(F_\infty,\ S \le c) + P(S > c)$$

$$\Delta P_L(c) = P(S>c)\,\big[1 - P(F_\infty \mid S>c)\big]$$

- $\Delta P_L(c)$ is the excess error caused by the cap.
- $P(F_\infty \mid S>c)$ is the fraction of slow shots that would have failed anyway. It is an interesting number in its own right: if it is close to 1, capping is nearly free.
- Without a fallback, an aborted decode counts as a failure: $P(F_c \mid S>c) = 1$.

**Requirement.** The full budget is $P_L(c) + P_{\text{ovf}} \le \varepsilon$, where $P_{\text{ovf}}$ is the buffer overflow probability. Note that $P(S>c)\,P(F_c\mid S>c) \le \varepsilon$ is necessary but not sufficient.

### Define the cap as a hardware budget, not a quantile
A cap set at the $(1-10^{-3})$ quantile produces a floor near $10^{-3}$ by definition, and a reviewer will say so. Instead, fix $c$ as an absolute operation budget (e.g. ops per µs at an assumed clock) and sweep $d$ and $p$.

For Relay-BP, the deployed small iteration and leg limits are one natural value of $c$, so the shipped configuration appears as a point on the curve.

**Probed caps.** Use three caps, $c_3, c_4, c_5$. Calibrate each once at a reference point $(d_{\text{ref}}, p_{\text{ref}})$ so that $P(S > c_j) \approx 10^{-j}$ there, i.e. cutoffs at roughly the $10^{-3}$, $10^{-4}$ and $10^{-5}$ tail levels. Then hold each $c_j$ fixed in absolute operations while sweeping $d$ and $p$, so the floor moves with the code instead of being fixed by definition.

### Key figures
1. **LER vs $d$ at fixed $c$** ($c_3, c_4, c_5$; fixed $p$). As $d$ grows, $p_L^{(0)}$ falls but window size and $P(S>c)$ rise, so the curve reaches a minimum and turns up. This optimal distance set by the decoder's tail is the headline figure.
2. **$P_L(c)$ vs $c$** at fixed $(d,p)$, with $p_L^{(0)}$ as an asymptote.
3. **$P(F_\infty \mid S>c)$ vs $c$**: are slow shots the failing shots?

### Budget
At most $10^6$ shots per $(d,p)$ point, for every setting. Only floors and cutoffs at $10^{-3}$, $10^{-4}$ and $10^{-5}$ are probed. At $10^{-4}$ that leaves about 100 capped shots per point (about 10% relative error); $c_5$ is marginal (about 10 shots at the reference point). No importance sampling is needed for this experiment. Data collection is batch-parallel: shots are split into fixed-size chunks with per-chunk seeds, decoded by $N$ worker processes (`scripts/pipeline.sh N`).

## 5. Experiment 2: Tail-class characterization

**Goal.** Determine which tail class each decoder's work distribution follows (light, lognormal, or Pareto), extract its parameters, and convert them into buffer-size scaling.

**Data.** Reuse the uncapped Experiment 1 data. No new decoding is required except possibly more shots at a few points.

### Candidate classes and what they imply for the buffer

**Spare servers.** With $k$ decoder servers and load $\rho = \lambda E[S]$ (measured in server units, where $\lambda$ is the window arrival rate), define the number of spare servers as

$$m \approx \lceil k - \rho \rceil,$$

the number of servers that can be stuck on slow decodes at once before the remaining servers can no longer keep up with arrivals. With heavy tails, overflow then requires about $m$ simultaneous slow decodes, so the per-window tail probability enters raised to the power $m$. For sliding-window decoding, $k=1$ and the D/G/1 column applies.

To leading order (single-big-jump asymptotics for D/G/1; Foss–Korshunov multi-jump results for D/G/k with out-of-order release):

| Tail class of $S$ | $P(S>x)$ | Sliding window (D/G/1) buffer | PWD (D/G/k) buffer |
|---|---|---|---|
| Light (exponential) | $e^{-\theta x}$ | $\sim \log(1/\varepsilon)$ | $\sim \log(1/\varepsilon)$, smaller constant |
| Lognormal $(\mu,\sigma)$ | $\bar\Phi\big((\ln x-\mu)/\sigma\big)$ | $e^{\mu+\sigma\sqrt{2\ln(1/\varepsilon)}}$ | $e^{\mu+\sigma\sqrt{2\ln(1/\varepsilon)/m}}$ |
| Pareto (index $\alpha>1$) | $(x/x_0)^{-\alpha}$ | $\sim \varepsilon^{-1/(\alpha-1)}$ | $\sim \varepsilon^{-1/(m(\alpha-1))}$ |

Here $\bar\Phi$ is the standard normal survival function, and $\theta, \mu, \sigma, x_0, \alpha$ are the fitted parameters of each class.

Takeaways:
- **The load barely matters.** The factor $\rho/(1-\rho)$ shifts $\ln x$ only slightly. Overflow comes from one huge decode, not gradual buildup, so overprovisioning throughput does not fix a heavy tail.
- **Lognormal at $\varepsilon=10^{-12}$:** $z\approx 7.0$, so the buffer is about $e^{7\sigma}$ median decodes. That is roughly 8× for $\sigma=0.3$, 33× for $\sigma=0.5$, and 1100× for $\sigma=1$. With $m=3$ spare servers, $\sigma=1$ drops to about 40×.
- **PWD changes the constant, not the class**, except that for a Pareto tail it divides the polynomial exponent by $m$. A fat tail is where parallel windows help most.

### Plots
For each $(d,p)$, plot the empirical complementary CDF $\bar F(x) = P(S>x)$ on log-log axes, without fitted curves.

### Fitting and tests
- **Peaks over threshold.** Fit a generalized Pareto distribution GPD$(\xi,\sigma)$ by maximum likelihood to the excesses above a high threshold $u$, and report the shape $\xi$ as the tail parameter:
  - $\xi>0$: a heavy tail, $\bar F(x)\sim x^{-1/\xi}$;
  - $\xi=0$: exponential;
  - $\xi<0$: bounded.

  Lognormal and Pareto are deliberately not distinguished: at any finite threshold both give $\xi>0$, and for buffer sizing the power-law reading is the conservative one.
- **Tests.** Only likelihood-ratio procedures:
  - a test of $\xi=0$ (exponential tail) against free $\xi$, with a $\chi^2_1$ null distribution;
  - a profile-likelihood (inverted likelihood-ratio) 95% confidence interval for $\xi$.

  Refit at a 10× higher threshold to check stability.
- **Censoring.** If any decodes hit the (relaxed) Relay-BP limit, use a censored likelihood. Better still, raise the limit until the censored fraction is far below the tail probabilities being fitted.

### Hypotheses
- **MWPM:** a lognormal fit has been reported (and observed by us) on the bulk; whether it holds in the far tail is the question.
- **Relay-BP:** BP non-convergence behaves like escape from traps, with escape time roughly geometric given trap depth. Mixing over trap depths can produce a **power-law (Pareto)** tail. Relay's leg structure may lighten it if legs act as nearly independent restarts. Evidence of long tails in BP-family iteration counts: Wang–Li–Mueller (HPCA'26, BP-SF), Landa (IBM), and the beam-search decoder paper (arXiv:2512.07057), which notes BP-LSD tails about 100× the mean.

### Outputs
1. Tail class and fitted parameters per (decoder, code, $d$, $p$), and how they trend with $d$ and $p$.
2. Required buffer size vs $\varepsilon$ (from $10^{-6}$ to $10^{-15}$) for sliding windows and for PWD with $m = 1, 2, 3$, using the table above with the measured parameters.

## 6. Conclusions the letter is set up to draw

- Capping at any practical work budget produces an LER floor that grows with $d$, so uncapped decoding plus a buffer is required.
- **MWPM:** if lognormal, the sliding-window buffer scales as $e^{\sigma\sqrt{2\ln(1/\varepsilon)}}$, and PWD with out-of-order release shrinks it by $\sqrt{m}$ in the exponent. That is likely necessary to keep the buffer manageable.
- **Relay-BP:**
  - light-tailed → buffer $\sim \log(1/\varepsilon)$, a significant positive result;
  - lognormal → same conclusion as MWPM;
  - Pareto → buffer polynomial in $1/\varepsilon$, and parallel windows for BP (not yet demonstrated for BB codes) become essential.

## 7. Suggested letter outline

1. Problem: syndrome streams, variable decode work, why a buffer is needed.
2. Capping creates an error floor (Experiment 1, headline figure).
3. Measured tail classes for MWPM and Relay-BP (Experiment 2).
4. Buffer scaling table, with buffer vs $\varepsilon$ for measured parameters.
5. Implications: PWD with out-of-order release, decoder switching (D/G/1 on the weak decoder), and the open question for BP.

## 8. Key references

- Skoric et al., *Parallel window decoding enables scalable fault tolerant quantum computation*, Nat. Commun. 14, 7040 (2023).
- Roffe et al., *Decoding across the quantum LDPC code landscape*, arXiv:2005.07016.
- Toshio et al., *Decoder switching: breaking the speed–accuracy tradeoff in real-time QEC*, arXiv:2510.25222.
- Müller et al., *Improved belief propagation is sufficient for real-time decoding of quantum memory* (Relay-BP), arXiv:2506.01779.
- Bravyi et al., *High-threshold and low-overhead fault-tolerant quantum memory*, Nature 627, 778 (2024).
- Gong, Cammerer, Renes, *Toward low-latency iterative decoding of QLDPC codes under circuit-level noise*, arXiv:2403.18901.
- Wang, Li, Mueller, *Fully parallelized BP decoding for quantum LDPC codes can outperform BP-OSD* (HPCA'26), arXiv:2507.00254.
- Landa, *Low-weight quantum syndrome errors in belief propagation decoding*, arXiv:2603.19126.
- Higgott & Gidney, *Sparse Blossom*, arXiv:2303.15933.
- Foss & Korshunov, heavy-tailed asymptotics for multi-server (GI/GI/k) queues; Veraverbeke, subexponential waiting-time asymptotics for GI/G/1.
- Clauset, Shalizi, Newman, *Power-law distributions in empirical data*, SIAM Review (2009).
